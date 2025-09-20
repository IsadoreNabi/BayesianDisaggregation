#' Compute likelihood vector from a prior matrix via SVD (center-only, robust)
#'
#' Builds a sectoral likelihood \eqn{L} (length \eqn{K}) from the prior weights
#' matrix \eqn{P \in \mathbb{R}^{T \times K}} by taking the absolute value of the
#' first right singular vector of the centered matrix (no scaling), then
#' normalizing to the unit simplex. Includes input validation, optional row
#' renormalization, and a safe fallback when PC1 is degenerate.
#'
#' @param P Numeric matrix \eqn{T \times K}; prior weights per period.
#' @param renormalize_rows Logical; if \code{TRUE} (default), rows of \code{P} that
#'   are within tolerance of summing to 1 are renormalized. Otherwise an error is thrown.
#' @param tol Numeric tolerance for simplex checks (default \code{1e-12}).
#'
#' @details
#' \strong{Validation:} \code{P} must be a finite, non-negative numeric matrix.
#' Each row must either (i) already sum to 1 within \code{tol} or (ii) be renormalizable
#' within \code{tol}. Rows with (near-)zero sums are not renormalizable and raise an error.
#' Missing values are not allowed.
#'
#' \strong{Algorithm (exactly as implemented):}
#' \enumerate{
#'   \item Center columns over time: \code{X <- scale(P, center = TRUE, scale = FALSE)}.
#'   \item Compute SVD: \code{sv <- svd(X)}.
#'   \item Take the first right singular vector (first column of V matrix); set \eqn{l = |v|}.
#'   \item If \eqn{\sum l \leq tol} or PC1 is degenerate, fall back to column means of \code{P} (over time) and renormalize.
#'   \item Otherwise, \eqn{L = l / \sum l}.
#' }
#'
#' @return Numeric vector \eqn{L} of length \eqn{K} (non-negative, sums to 1). Attributes:
#' \itemize{
#'   \item \code{"pc1_loadings"}: signed PC1 loadings \eqn{v}.
#'   \item \code{"explained_var"}: fraction of variance explained by PC1.
#'   \item \code{"fallback"}: \code{TRUE} if column-mean fallback was used.
#' }

#'
#' @examples
#' set.seed(123)
#' T <- 10; K <- 4
#' P <- matrix(rexp(T*K), nrow = T); P <- P / rowSums(P)
#' L <- compute_L_from_P(P)
#' stopifnot(length(L) == K, all(L >= 0), abs(sum(L) - 1) < 1e-12)
#'
#' @seealso \code{\link{spread_likelihood}}, \code{\link{bayesian_disaggregate}}
#' @export
compute_L_from_P <- function(P, renormalize_rows = TRUE, tol = 1e-12) {
  # --- validation ---
  if (!is.matrix(P) || !is.numeric(P)) stop("P must be a numeric matrix.")
  if (any(!is.finite(P))) stop("P contains non-finite values.")
  if (any(P < 0)) stop("P contains negative values; expected non-negative simplex rows.")

  # quick path K = 1
  if (ncol(P) == 1L) {
    L <- 1
    attr(L, "pc1_loadings") <- 1
    attr(L, "explained_var") <- 1
    attr(L, "fallback") <- FALSE
    return(L)
  }

  rs <- rowSums(P)
  if (any(rs <= tol)) {
    stop("At least one row of P has (near-)zero sum; cannot renormalize.")
  }
  if (any(abs(rs - 1) > tol)) {
    if (renormalize_rows) {
      # row-wise normalization
      P <- sweep(P, 1, rs, "/")
    } else {
      stop("Row sums of P are not 1 within 'tol' and renormalize_rows = FALSE.")
    }
  }

  # --- center only, then SVD ---
  X  <- scale(P, center = TRUE, scale = FALSE)
  sv <- svd(X)

  d2 <- sv$d^2
  explained_var <- if (sum(d2) > 0) d2[1] / sum(d2) else 0

  v <- sv$v[, 1]
  l <- abs(v)
  fallback <- FALSE

  if (sum(l) <= tol || explained_var <= tol || any(!is.finite(l))) {
    # fallback: column means over time, then normalize
    l <- colMeans(P)
    if (any(!is.finite(l)) || sum(l) <= tol) {
      stop("Degenerate PC1 and invalid column means; cannot compute L.")
    }
    fallback <- TRUE
  }

  L <- l / sum(l)
  attr(L, "pc1_loadings") <- v
  attr(L, "explained_var") <- explained_var
  attr(L, "fallback") <- fallback
  L
}

#' Spread a likelihood vector across time with a chosen temporal pattern
#'
#' Expands a sectoral likelihood \eqn{L} (length \eqn{K}) into a \eqn{T \times K}
#' matrix by applying a temporal weight profile and then row-normalizing to the
#' simplex.
#'
#' @param L Numeric vector (length \eqn{K}); sectoral likelihood. It is normalized
#'   internally to sum to 1 before spreading.
#' @param T_periods Integer; number of time periods \eqn{T}.
#' @param pattern Temporal pattern for the weights across time; one of
#'   \code{"constant"}, \code{"recent"}, \code{"linear"}, \code{"bell"}.
#'
#' @details
#' Given \eqn{L} and a non-negative time-weight vector \eqn{w_t}, the function
#' replicates \eqn{L} across rows and applies \eqn{w} elementwise, then
#' \emph{row-normalizes} using \code{row_norm1()}:
#' \deqn{LT_{t,k} \propto w_t \cdot L_k, \qquad \sum_k LT_{t,k} = 1 \ \forall t.}
#'
#' Patterns:
#' \itemize{
#'   \item \code{"constant"}: \eqn{w_t = 1}.
#'   \item \code{"recent"}: linearly increasing in \eqn{t} (more weight to later periods).
#'   \item \code{"linear"}: affine ramp from first to last period.
#'   \item \code{"bell"}: symmetric bell-shaped profile centered at \eqn{T/2}.
#' }
#'
#' @return Numeric matrix \eqn{T \times K}; each row sums to 1.
#'
#' @examples
#' set.seed(1)
#' K <- 5; T <- 8
#' L <- runif(K); L <- L / sum(L)
#' LT <- spread_likelihood(L, T, "recent")
#' stopifnot(nrow(LT) == T, ncol(LT) == K, all(abs(rowSums(LT) - 1) < 1e-12))
#'
#' @seealso \code{\link{compute_L_from_P}}, \code{\link{bayesian_disaggregate}}
#' @export
spread_likelihood <- function(L, T_periods,
                              pattern = c("constant", "recent", "linear", "bell")) {
  # --- validation ---
  if (!is.numeric(L) || length(L) < 1L || any(!is.finite(L))) {
    stop("L must be a finite numeric vector of length >= 1.")
  }
  if (any(L < 0)) stop("L must be non-negative.")
  if (length(T_periods) != 1L || !is.numeric(T_periods) || T_periods < 1) {
    stop("T_periods must be a positive integer.")
  }
  T_periods <- as.integer(T_periods)

  pattern <- match.arg(pattern)

  # normalize L to simplex
  sL <- sum(L)
  if (sL <= 0) stop("Sum(L) must be > 0.")
  L <- as.numeric(L / sL)

  # time weights
  t <- seq_len(T_periods)
  w <- switch(
    pattern,
    constant = rep(1, T_periods),
    recent   = 0.5 + 0.5 * (t / T_periods),
    linear   = 0.3 + 1.4 * ((t - 1) / max(1, T_periods - 1)),
    bell     = exp(- (t - T_periods/2)^2 / (2 * (T_periods/4)^2))
  )

  if (any(!is.finite(w)) || any(w < 0)) stop("Temporal weights became invalid.")
  if (all(w == 0)) stop("Temporal weights are all zero.")

  # replicate and apply weights, then row-normalize
  LT <- matrix(rep(L, each = T_periods), nrow = T_periods) * w
  row_norm1(LT)
}