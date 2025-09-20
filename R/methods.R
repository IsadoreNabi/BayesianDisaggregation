#' Weighted-average posterior (convex combination)
#'
#' Computes \eqn{W = \lambda P + (1-\lambda)\,LT} row-wise and renormalizes each
#' row to the unit simplex (suma 1). Expects matrices no negativas con las mismas
#' dimensiones \eqn{T \times K}.
#'
#' @param P Prior matrix (\eqn{T \times K}); filas no negativas.
#' @param LT Likelihood matrix (\eqn{T \times K}); filas no negativas.
#' @param lambda Mixing weight in \eqn{[0,1]} (default \code{0.7}).
#'
#' @return Posterior matrix \eqn{W} (\eqn{T \times K}), filas suman 1.
#'
#' @examples
#' T <- 5; K <- 3
#' P  <- matrix(runif(T*K), T); P  <- P / rowSums(P)
#' LT <- matrix(runif(T*K), T); LT <- LT / rowSums(LT)
#' W <- posterior_weighted(P, LT, 0.6)
#' stopifnot(all(abs(rowSums(W)-1) < 1e-12))
#'
#' @seealso \code{\link{posterior_multiplicative}}, \code{\link{posterior_dirichlet}},
#'   \code{\link{posterior_adaptive}}
#' @export
posterior_weighted <- function(P, LT, lambda = 0.7) {
  if (!is.matrix(P) || !is.matrix(LT)) stop("P and LT must be matrices.")
  if (!all(dim(P) == dim(LT))) stop("P and LT must have the same dimensions.")
  if (!is.numeric(lambda) || lambda < 0 || lambda > 1) stop("lambda must be in [0,1].")
  if (any(P < 0) || any(LT < 0)) stop("P and LT must be non-negative.")
  row_norm1(lambda * P + (1 - lambda) * LT)
}

#' Multiplicative posterior (Hadamard product + renormalization)
#'
#' Computes \eqn{W \propto P \odot LT}{W ~ P * LT} (producto elemento a elemento) y
#' renormaliza cada fila a la simplex. Útil cuando prior y likelihood deben
#' reforzarse mutuamente.
#'
#' @param P Prior matrix (\eqn{T \times K}); filas no negativas.
#' @param LT Likelihood matrix (\eqn{T \times K}); filas no negativas.
#'
#' @return Posterior matrix \eqn{W} (\eqn{T \times K}), filas suman 1.
#'
#' @examples
#' T <- 4; K <- 4
#' P  <- matrix(runif(T*K), T); P  <- P / rowSums(P)
#' LT <- matrix(runif(T*K), T); LT <- LT / rowSums(LT)
#' W <- posterior_multiplicative(P, LT)
#'
#' @seealso \code{\link{posterior_weighted}}, \code{\link{posterior_dirichlet}},
#'   \code{\link{posterior_adaptive}}
#' @export
posterior_multiplicative <- function(P, LT) {
  if (!is.matrix(P) || !is.matrix(LT)) stop("P and LT must be matrices.")
  if (!all(dim(P) == dim(LT))) stop("P and LT must have the same dimensions.")
  if (any(P < 0) || any(LT < 0)) stop("P and LT must be non-negative.")
  row_norm1(P * LT)
}

#' Dirichlet-conjugate posterior (analytical mean)
#'
#' Interpreta \eqn{P} y \eqn{LT} como proporciones convertidas a pseudo-cuentas
#' con concentración total \eqn{\alpha_{base}=1/\gamma}. El posterior medio es
#' \eqn{W = (P\,\alpha_{base} + LT\,\alpha_{base}) / \mathrm{rowSums}(\cdot)}.
#'
#' @param P Prior matrix (\eqn{T \times K}); filas no negativas que suman 1 (o cercanas).
#' @param LT Likelihood matrix (\eqn{T \times K}); filas no negativas que suman 1 (o cercanas).
#' @param gamma Positive uncertainty factor (default \code{0.1}); menor \eqn{\gamma}
#'   implica mayor concentración (más “seguro”).
#'
#' @return Posterior matrix \eqn{W} (\eqn{T \times K}), filas suman 1.
#'
#' @examples
#' T <- 6; K <- 3
#' P  <- matrix(runif(T*K), T);  P  <- P  / rowSums(P)
#' LT <- matrix(runif(T*K), T); LT <- LT / rowSums(LT)
#' W <- posterior_dirichlet(P, LT, gamma = 0.1)
#'
#' @seealso \code{\link{posterior_weighted}}, \code{\link{posterior_multiplicative}},
#'   \code{\link{posterior_adaptive}}
#' @export
posterior_dirichlet <- function(P, LT, gamma = 0.1) {
  if (!is.matrix(P) || !is.matrix(LT)) stop("P and LT must be matrices.")
  if (!all(dim(P) == dim(LT))) stop("P and LT must have the same dimensions.")
  if (any(P < 0) || any(LT < 0)) stop("P and LT must be non-negative.")
  if (!is.numeric(gamma) || gamma <= 0) stop("gamma must be > 0.")
  alpha_base <- 1 / gamma
  A <- P * alpha_base + LT * alpha_base
  A / rowSums(A)
}

#' Adaptive posterior based on sector volatility
#'
#' Ajusta el *mixing* por sector según la volatilidad histórica del prior:
#' \deqn{\phi_k = \min\!\left(\frac{\sigma_k/\mu_k}{\overline{\sigma/\mu}}, 0.8\right)},
#' donde \eqn{\sigma_k/\mu_k} es la desviación estándar relativa de la columna \eqn{k}.
#' Sectores más volátiles reciben más peso de \eqn{LT}. Devuelve
#' \eqn{W = (1-\phi)\odot P + \phi\odot LT}{W = (1-phi)*P + phi*LT} renormalizado por filas.
#'
#' @param P Prior matrix (\eqn{T \times K}); filas no negativas.
#' @param LT Likelihood matrix (\eqn{T \times K}); filas no negativas.
#'
#' @return Posterior matrix \eqn{W} (\eqn{T \times K}), filas suman 1.
#'
#' @examples
#' set.seed(1)
#' T <- 8; K <- 5
#' P  <- matrix(runif(T*K), T);  P  <- P  / rowSums(P)
#' LT <- matrix(runif(T*K), T); LT <- LT / rowSums(LT)
#' W <- posterior_adaptive(P, LT)
#'
#' @seealso \code{\link{posterior_weighted}}, \code{\link{posterior_multiplicative}},
#'   \code{\link{posterior_dirichlet}}
#' @export
posterior_adaptive <- function(P, LT) {
  if (!is.matrix(P) || !is.matrix(LT)) stop("P and LT must be matrices.")
  if (!all(dim(P) == dim(LT))) stop("P and LT must have the same dimensions.")
  if (any(P < 0) || any(LT < 0)) stop("P and LT must be non-negative.")

  # sd/mean robustificado: si mean==0 => 0
  rel_sd <- function(x) {
    m <- mean(x)
    if (!is.finite(m) || m <= 0) return(0)
    s <- sd(x)
    if (!is.finite(s)) return(0)
    s / m
  }

  sigma_k <- apply(P, 2, rel_sd)
  bar_sigma <- mean(sigma_k)
  if (!is.finite(bar_sigma) || bar_sigma <= 0) {
    # Sin señal de volatilidad: usa weighted con lambda=1 (i.e., prior)
    return(row_norm1(P))
  }

  phi_k <- pmin(sigma_k / bar_sigma, 0.8)
  row_norm1(sweep(P,  2, (1 - phi_k), `*`) +
            sweep(LT, 2,      phi_k , `*`))
}