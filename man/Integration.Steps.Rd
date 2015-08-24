\encoding{latin1}
\name{Integration.Steps}
\alias{Integration.Steps}
\title{Subintervals and Integration Steps To Approximate First-Passage-Time Densities}
\description{
  According to the First-Passage-Time Location (FPTL) function and the arguments in the function call, this function  
  calculates suitable subintervals and integration steps in order to approximate the first-passage-time (f.p.t.) density. 
}
\usage{
Integration.Steps(sfptl, variableStep = TRUE, from.t0 = FALSE,
                  to.T = FALSE, n = 250, p = 0.2, alpha = 1)
}
\arguments{
  \item{sfptl}{an object of class \dQuote{summary.fptl}.}
  \item{variableStep}{a logical value indicating whether a variable integration step is used.}
  \item{from.t0}{a logical value indicating whether the approximation should be calculated from the lower end of the
interval considered, \eqn{t_0}{t0}, specified in the object used as \code{sfptl} argument.}
  \item{to.T}{a logical value indicating whether the approximation should be calculated to the upper end of the
interval considered, \eqn{T}, specified in the object used as \code{sfptl} argument.}
	\item{n}{Number of points used to determine the integration step in subintervals \eqn{[t_i^*, t_{max,i}^+]}{[t[i]*, tmax[i]^+]}, 
\eqn{i=1, \ldots, m}{i=1,..., m}, from interesting instants provided by the FPTL function.}
  \item{p}{Ratio of n used to determine the integration step in subintervals \eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]}, 
  \eqn{i=1, \ldots, m}{i=1,..., m}, \eqn{[t_0, t_1^*]}{[t0, t[1]*]} and \eqn{[t_{max,m}^{+}, T]}{[tmax[m]^+, T]}.}
  \item{alpha}{Parameter used to determine the integration step in subintervals 
\eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]}, \eqn{i=1, \ldots, m}{i=1,..., m}, \eqn{[t_0, t_1^*]}{[t0, t[1]*]} and 
\eqn{[t_{max,m}^{+}, T]}{[tmax[m]^+, T]}, in order to reduce the computational cost of approximating the f.p.t. density function 
in those cases where \eqn{t_{i+1}^* - t_{max,i}^+ >> t_{max,i}^{-} - \thinspace t_i^*}{t[i+1]* - tmax[i]^+ >> tmax[i]^- - t[i]*}, 
for some \eqn{i}, \eqn{t_1^* - t_0 >> t_{max,1}^{-} - \thinspace t_1^*}{t[1]* - t0 >> tmax[1]^- - t[1]*} or 
\eqn{T - t_{max,m}^+ >> t_{max,m}^{-} - \thinspace t_m^*}{T - tmax[m]^+ >> tmax[m]^- - t[m]*}, respectively.}
} 
\details{
Based on the information provided by the FPTL function contained in the \code{sfptl} object, this function computes and returns 
suitable subintervals and integration steps in order to approximate the density 
function of the f.p.t. variable according to the other arguments in the function call. \cr

When the \code{sfptl} object is of length greater than 1, it comes from an unconditioned f.p.t. problem. Each component is associated with 
the same f.p.t. problem conditioned on different values of the initial distribution
(equally spaced in the range of the distribution). Let \eqn{x_{0,j}}{x0[j]}, \eqn{j=1, \ldots, N}, these values. For each initial value
\eqn{x_{0,j}}{x0[j]} let \eqn{t_{i,j}^*}{t[i, j]^*}, \eqn{t_{max,i,j}^-}{tmax[i, j]^-} and \eqn{t_{max,i,j}^+}{tmax[i, j]^+}, \eqn{i=1,\ldots \ m_j}{i=1,...,m[j]}, 
the interesting time instants provided by the FPTL function and stored in the \code{instants} component of the j-th list in the \code{sfptl} object.
Then, the time instants \eqn{\{t_{i,j}, i=1, 2, \ldots, 2m_j \}}{{t[i, j], i=1,2,...,2m[j]}}, where 
\ifelse{latex}{\deqn{t_{i,j} = \left\{ \begin{array}{ll} t_{(i+1)/2, \thinspace j}^* & \mbox{for i odd} \smallskip \\ t_{max, \thinspace i/2, \thinspace j}^+ & \mbox{for i even} \end{array} 
\right. ,}}{\deqn{t[i, j] = t[(i+1)/2, j]*,} for i odd, and \deqn{t[i, j] = tmax[i/2, j]^+,} for i even,}
provide a suitable partition of interval \eqn{[t_0, T]}{[t0, T]} to approximate the f.p.t density
for the fixed value \eqn{x_{0,j}}{x0[j]} of the initial distribution. \cr

If the \code{sfptl} object is of length 1, it comes from a conditioned f.p.t. problem. In this case we denote the interesting time 
instants provided by the FPTL function and stored in the \code{sfptl} object by \eqn{t_{i,1}^*}{t[i, 1]^*}, \eqn{t_{max,i,1}^-}{tmax[i, 1]^-} 
and \eqn{t_{max,i,1}^+}{tmax[i, 1]^+}. \cr

In what follows, \eqn{\lceil x \rceil}{[x]} is the \ifelse{latex}{smallest integer not less than \eqn{x}}{integer part of \eqn{x}}.

For each list in the \code{sfptl} object the function computes
\deqn{h_{i,j} = \frac{t_{max,i,j}^{+} - \thinspace t_{i,j}^*}{n_{i,j}} , i=1, \ldots, m_j,}{h[i,j] = (tmax[i, j]^+ - t[i, j]*)/n[i, j], i=1,...,m[j],} 
where
\deqn{n_{i,j} = \lceil n k_{i,j} \rceil}{n[i, j] = [n k[i, j]] + 1}
and
\deqn{k_{i,j} = \frac{t_{max,i,j}^+ - t_{i,j}^{*}}{t_{max,i,j}^{-} - \thinspace t_{i,j}^*} \ .}{k[i, j] = (tmax[i, j]^+ - t[i, j]*)/(tmax[i, j]^- - t[i, j]*).} \cr 

If \code{variableStep = TRUE} also computes, for
\if{latex}{\eqn{0 < p, \alpha \leq 1}:}
\if{html}{\emph{0 < p, alpha \out{&le;} 1}:}
\if{text}{\out{0 < p, alpha <= 1}:}
\itemize{
\item \eqn{h_{i,j}^* = \displaystyle{\frac{t_{i+1,j}^* - t_{max,i,j}^+}{n_{i,j}^*}} , i=1, \ldots, m_j-1,}{h[i, j]* = (t[i+1, j]* - tmax[i, j]^+)/n[i, j]*, i=1,...,m[j]-1,}
where
\ifelse{latex}{\deqn{n_{i,j}^* = \left\{ \begin{array}{ccc}
          \lceil n \thinspace p \thinspace k_{i,j}^* \rceil &  \mbox{if} & k_{i,j}^* \leq 1 \vspace{0.25cm} \\
          \lceil n \thinspace p \thinspace k_{i,j}^{*^{alpha}} \rceil & \mbox{if} & k_{i,j}^* > 1
        \end{array}
\right.} and}{\deqn{n[i, j]* = [n p k[i, j]*] + 1 ,} if \eqn{k[i, j]* \le 1}, and \deqn{n[i, j]* = [n p k[i, j]*^alpha] + 1 ,} if \eqn{k[i, j]* > 1}, with}
\deqn{k_{i,j}^* = \frac{t_{i+1,j}^* - \thinspace t_{max,i,j}^+}{t_{max,i,j}^{-} - \thinspace t_{i,j}^*} \ .}{k[i, j]* = (t[i+1, j]* - tmax[i, j]^+)/(tmax[i, j]^- - t[i, j]*) .}

If \eqn{h_{i,j}^* < h_{i,j}}{h[i, j]* < h[i, j]}, we then set \eqn{t_{max,i,j}^+}{tmax[i, j]^+} equal to \eqn{t_{i+1,j}^*}{t[i+1, j]*} and \eqn{h_{i,j}}{h[i, j]} is recalculated.
 
\item \eqn{h_{0,j}^* = \displaystyle{\frac{t_{1,j}^* - t_0}{n_{0,j}^*}} , }{h[0, j]* = (t[1, j]* - t0)/n[0, j]* ,}
where
\ifelse{latex}{\deqn{n_{0,j}^* = \left\{ \begin{array}{ccc}
          \lceil n \thinspace p \thinspace k_{0,j}^* \rceil &  \mbox{if} & k_{0,j}^* \leq 1 \vspace{0.25cm} \\
          \lceil n \thinspace p \thinspace k_{0,j}^{*^{alpha}} \rceil & \mbox{if} & k_{0,j}^* > 1
        \end{array}
 \right.} and}{\deqn{n[0, j]* = [n p k[0, j]*] + 1 ,} if \eqn{k[0, j]* \le 1}, and \deqn{n[0, j]* = [n p k[0, j]*^alpha] + 1,} if \eqn{k[0, j]* > 1}, with}
\deqn{k_{0,j}^* = \frac{t_{1,j}^* - \thinspace t_0}{t_{max,1,j}^{-} - \thinspace t_{1,j}^*} \ ,}{k[0, j]* = (t[1, j]* - t0)/(tmax[1, j]^- - t[1, j]*) ,} 

when the \code{sfptl} object is of length 1 and \code{from.t0 = TRUE}, or the \code{sfptl} object is of length greater than 1.

If \eqn{h_{0,j}^* < h_{1,j}}{h[0, j]* < h[1, j]}, we then set \eqn{t_{1,j}^*}{t[1, j]*} equal to  \eqn{t_0}{t0}  and \eqn{h_{1,j}}{h[1, j]} is recalculated.

\item \eqn{h_{m_j,j}^* = \displaystyle{\frac{T - \thinspace t_{max,m_j,j}^+}{n_{m_j,j}^*}} , }{h[m[j], j]* = (T - tmax[m[j], j]^+)/n[m[j], j]* ,}
where
\ifelse{latex}{\deqn{n_{m_j,j}^* = \left\{ \begin{array}{ccc}
          \lceil n \thinspace p \thinspace k_{m_j,j}^* \rceil &  \mbox{if} & k_{m_j,j}^* \leq 1 \vspace{0.25cm} \\
          \lceil n \thinspace p \thinspace k_{m_j,j}^{*^{alpha}} \rceil & \mbox{if} & k_{m_j,j}^* > 1
        \end{array}
 \right. } and}{\deqn{n[m[j], j]* = [n p k[m[j], j]*] + 1 ,} if \eqn{k[m[j], j]* \le 1}, and \deqn{n[m[j], j]* = [n p k[m[j], j]*^alpha] + 1 ,} if \eqn{k[m[j], j]* > 1}, with}
\deqn{k_{m_j,j}^* = \frac{T - \thinspace t_{max,m_j,j}^+}{t_{max,m_j,j}^{-} - \thinspace t_{m_j,j}^*} \ ,}{k[m[j], j]* = (T - tmax[m[j], j]^+)/(tmax[m[j], j]^- - t[m[j], j]^*) ,} 
when the \code{sfptl} object is of length 1 and \code{to.T = TRUE}, or the \code{sfptl} object is of length greater than 1.

If \eqn{h_{m_j,j}^* < h_{m_j,j}}{h[m[j], j]* < h[m[j], j]}, we then set \eqn{t_{max,m_j,j}^+}{tmax[m[j], j]^+} equal to \eqn{T} and \eqn{h_{m_j,j}}{h[m[j], j]} is recalculated. \cr
}

\if{latex}{\eqn{p \geq 0.1} and \eqn{0.75 \leq alpha \leq 1} are recommended; otherwise, some integration steps can be excessively large. \cr}
\if{html}{\emph{p} \out{&ge;} 0.1 and 0.75 \out{&le;} \emph{alpha} \out{&le;} 1 are recommended; if not, some integration steps can be excessively large. \cr} 
\if{text}{\out{p >= 0.1 and 0.75 <= alpha <= 1 are recommended; if not, some integration steps can be excessively large.} \cr} 

If the \code{sfptl} object is of length 1 (conditioned f.p.t. problem), the suitable subintervals and integration steps that the function provides are:
\itemize{
\item If \code{variableStep = TRUE},
\itemize{
\item \ifelse{latex}{\eqn{h_{i,1}} in subintervals \eqn{[t_{i,1}^*, t_{max,i,1}^+]}, \eqn{i=1, \ldots, m_1}.}{\emph{h[i, 1]} in subintervals 
\eqn{[t[i, 1]*, tmax[i, 1]^+]}, \eqn{i=1,...,m[1]}.}
\item \ifelse{latex}{\eqn{h_{i,1}^*} in subintervals \eqn{[t_{max,i,1}^+, t_{i+1,1}^*]}, \eqn{i=1, \ldots, m_1-1}. 
In these subintervals is possible to avoid applying the numerical algorithm to approximate the f.p.t. density provided that the
value of the approximate density at the time instant \eqn{t_{max,i,1}^+} is almost 0.}{\emph{h[i, 1]*} in subintervals \eqn{[tmax[i, 1]^+, t[i+1, 1]*]}, 
\eqn{i=1,...,m[1]-1}. In this subintervals is possible to avoid applying the numerical algorithm to approximate the f.p.t. density provided that the
value of the approximate density at the time instant \eqn{tmax[i,1]^+} is almost 0.}
\item \ifelse{latex}{\eqn{h_{0,1}^*} in subinterval \eqn{[t_0, t_{1,1}^*]}, if \code{from.t0 = TRUE}.}{\emph{h[0,1]*} in subinterval 
\eqn{[t0, t[1, 1]*]}, if \code{from.t0 = TRUE}.}
\item \ifelse{latex}{\eqn{h_{m_1,1}^*} in subinterval \eqn{[t_{max,m_1,1}^{+}, T]}, if \code{to.T = TRUE}.}{\emph{h[m[1], 1]*} in subinterval 
\eqn{[tmax[m[1], 1]^+, T]}, if \code{to.T = TRUE}.} \cr
}
\item If \code{variableStep = FALSE} the function computes 
\deqn{h = min \left\{ h_{i,1} \ , \ i=1, \ldots, m_1 \right\} \thinspace .}{h = min {h[i, 1] , i=1,...,m[1]}.}
Then
\itemize{
\item If \code{from.t0 = FALSE} and \code{to.T = FALSE}, \eqn{h} is readjusted to exactly split the interval \eqn{[t_{1,1}^*, t_{max,m_1,1}^+]}{[t[1, 1]*, tmax[m[1], 1]^+]}.
\item If \code{from.t0 = TRUE} and \code{to.T = FALSE}, \eqn{h} is readjusted to exactly split the interval \eqn{[t_0, t_{max,m_1,1}^+]}{[t0, tmax[m[1], 1]^+]}.
\item If \code{from.t0 = FALSE} and \code{to.T = TRUE}, \eqn{h} is readjusted to exactly split the interval \eqn{[t_{1,1}^*, T]}{[t[1, 1]*, T]}.
\item If \code{from.t0 = TRUE} and \code{to.T = TRUE}, \eqn{h} is readjusted to exactly split the interval \eqn{[t_0, T]}{[t0, T]}.
}

\emph{h} is a suitable fixed integration step in subintervals \eqn{[t_{i,1}^*, t_{max,i,1}^+]}{[t[i, 1]*, tmax[i, 1]^+]}, \eqn{i=1, \ldots, m_1}{i=1,...,m[1]}, in 
\eqn{[t_{max,i,1}^+, t_{i+1,1}^*]}{[tmax[i, 1]^+, t[i+1, 1]*]}, \eqn{i=1, \ldots, m_1-1}{i=1,...,m[1]-1}; 
in subintervals \eqn{[t_0, t_{1,1}^*]}{[t0, t[1, 1]*]} if \code{from.t0 = TRUE}, and in
\eqn{[t_{max,m_1,1}^{+}, T]}{[tmax[m[1], 1]^+, T]} if \code{to.T = TRUE}. The endpoints of such subintervals 
are readjusted according to this integration step. \cr
}

If the \code{sfptl} object is a list of length greater than 1 (unconditioned f.p.t problem), a common partition of 
the interval \eqn{[t_0, T]}{[t0, T]} is calculated from the suitable partitions of this interval for each 
fixed value of the initial distribution. \cr

Let, in unified form, \eqn{H_{r,j}}{H[r, j]}, \eqn{r=1, \ldots, 2m_j+1}{r=1,...,2m[j]+1}, the suitable integration steps
(calculated for each \eqn{j} in similar manner to the case of the \code{sfptl} object is of length 1) in 
subintervals \eqn{I_{r,j} = [t_{r-1,j}, t_{r,j}]}{I[r, j] = [t[r-1, j], t[r,j]]}, with \eqn{t_{0,j}=t_0}{t[0, j] = t0} and 
\eqn{t_{2m_j+1,j} = T}{t[2m[j]+1, j] = T}, \if{latex}{\eqn{\forall}} \if{html}{\out{&forall;}} \eqn{j=1, \ldots, N}{j=1,...,N}. 
Then, the ordered values of all time instants in the suitable partitions,
\eqn{t_{(1)}, \ldots, t_{(M)}}{t[(1)],...,t[(M)]}, provide a common suitable partition of the interval \eqn{[t_0, T]}{[t0, T]} 
in subintervals \eqn{[t_{(i-1)}, t_{(i)}]}{[t[(i-1)], t[(i)]]}, \eqn{i=1, \ldots, M \negthinspace + \negthinspace 1}{i=1,...,M+1},
where \eqn{t_{(0)} = t_0}{t[(0)] = t0} and \eqn{t_{(M \negthinspace + \negthinspace 1)} = T}{t[(M+1)] = T}. \cr

For this partition, the function computes
\itemize{ 
\item \if{latex}{\eqn{H_{i} = \mbox{min} \{ H_{r,j} : j=1, \ldots, N, \ and \ [t_{(i-1)}, t_{(i)}] \subseteq I_{r,j} \} , \ i=2, \ldots, M.}}
\if{html}{\emph{H[i] = min \{H[r, j] : j=1,...,N and [t[(i-1)], t[(i)]]\out{&sube;} I[r, j]\}, i=2,...,M.}}
\if{text}{H[i] = min {H[r, j] : j=1,...,N and [t[(i-1)], t[(i)]] is a subset of I[r, j]}, i=1,...,M.}
\item \if{latex}{\eqn{H_{1} = \mbox{min} \{ H_{r,j} : j=1, \ldots, N, \ and \ [t_0, t_{(1)}] \subseteq I_{r,j} \}}, if \code{from.t0 = TRUE}.}
\if{html}{\emph{H[1] = min \{H[r, j] : j=1,...,N and [t0, t[(1)]]\out{&sube;} I[r, j]\}}, if \code{from.t0 = TRUE}.}
\if{text}{H[1] = min {H[r, j] : j=1,...,N and [t0, t[(1)]] is a subset of I[r, j]}, if \code{from.t0 = TRUE}.}
\item \if{latex}{\eqn{H_{M \negthinspace + \negthinspace 1} = \mbox{min} \{ H_{r,j} : j=1, \ldots, N, \ and \ [t_{(M)}, T] \subseteq I_{r,j} \}}, 
if \code{to.T = TRUE}. \cr}
\if{html}{\emph{H[M+1] = min \{H[r, j] : j=1,...,N and [t[(M)], T]\out{&sube;} I[r, j]\}}, if \code{from.t0 = TRUE}.}
\if{text}{H[M+1] = min {H[r, j] : j=1,...,N and [t[(M)], T] is a subset of I[r, j]}, if \code{from.t0 = TRUE}.}
}

Thus, 
\itemize{
\item If \code{variableStep = TRUE}, the suitable subintervals and integrations steps that the function provides are
\itemize{
\item \ifelse{latex}{\eqn{H_i}}{\emph{H[i]}} in subintervals \eqn{[t_{(i-1)}, t_{(i)}]}{[t[(i-1)], t[(i)]]}, \eqn{i=2, \ldots, M}{i=2,...,M} 
\item \ifelse{latex}{\eqn{H_1}}{\emph{H[1]}} in subinterval \eqn{[t_0, t_{(1)}]}{[t0, t[(1)]]}, if \code{from.t0 = TRUE}.
\item \ifelse{latex}{\eqn{H_{M \negthinspace + \negthinspace 1}}}{\emph{H[M+1]}} in subinterval \eqn{[t_{(M)}, T]}{[t[(M)], T]}, if \code{to.T = TRUE}.
}
Each integration step is readjusted to exactly split the corresponding subinterval. \cr
\item If \code{variableStep = FALSE}, a suitable fixed integration step for any subinterval \eqn{[t_{(i-1)}, t_{(i)}]}{[t[(i-1)], t[(i)]]} is
\deqn{h = \mbox{min} \left\{ H_{i} : \begin{array}{l} i=1, \ldots, M \negthinspace + \negthinspace 1 \ \ \mbox{and it is not possible to avoid applying the} 
\\ \mbox{approximation algorithm in } [t_{(i-1)}, t_{(i)}] \ \forall \ x_{0,j}, j=1, \ldots, N. \end{array} \right\} . \ }{h = min {H[i] : i=1,...,M+1
and it is not possible to avoid applying the approximation algorithm in [t[(i-1)], t[(i)]] for all x0[j], j=1,...,N}}
Then
\itemize{
\item If \code{from.t0 = FALSE} and \code{to.T = FALSE}, \eqn{h} is readjusted to exactly split the interval \eqn{[t_{(1)}, t_{(M)}]}{[t[(1)], t[(M)]]}.
\item If \code{from.t0 = TRUE} and \code{to.T = FALSE}, \eqn{h} is readjusted to exactly split the interval \eqn{[t_0, t_{(M)}]}{[t0, t[(M)]]}.
\item If \code{from.t0 = FALSE} and \code{to.T = TRUE}, \eqn{h} is readjusted to exactly split the interval \eqn{[t_{(1)}, T]}{[t[(1)], T]}.
\item If \code{from.t0 = TRUE} and \code{to.T = TRUE}, \eqn{h} is readjusted to exactly split the interval \eqn{[t_0, T]}{[t0, T]}.
}

\emph{h} is a suitable fixed integration step in subintervals \eqn{[t_{(i-1)}, t_{(i)}]}{[t[(i-1)], t[(i)]]}, \eqn{i=2, \ldots, M}{i=2,..., M},  
in subintervals \eqn{[t_0, t_{(1)}]}{[t0, t[(1)]]} if \code{from.t0 = TRUE}, and in
\eqn{[t_{(M)}, T]}{[t[(M)], T]} if \code{to.T = TRUE}. The endpoints of such subintervals 
are readjusted according to this integration step. \cr
}
}
\value{
  A two-component list:
  \item{H}{A matrix of subintervals and integrations steps that we must consider in order to approximate the f.p.t. 
  density according to the information contained in the \code{sfptl} object and the arguments in the function call.}
  \item{skip}{A list of logical vectors indicating, for each subinterval, the values of the initial distribution for which 
  we must check whether it is possible to avoid applying the numerical algorithm.
  }  
}
\references{
Román, P., Serrano, J. J., Torres, F. (2008) First-passage-time location function: Application to determine first-passage-time densities in diffusion processes. \emph{Comput. Stat. Data Anal.}, \bold{52}, 4132--4146.

P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2012) An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.

P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2014) More general problems on first-passage times for diffusion processes: A new version of the fptdApprox R package. \emph{Applied Mathematics and Computation}, \bold{244}, 432--446.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
\code{\link{Approx.fpt.density}} to approximate f.p.t. densities from objects of class \dQuote{summary.fptl} and create objects of class \dQuote{fpt.density}.

\code{\link{summary.fptl}} to locate the f.p.t. variable and create objects of class \dQuote{summary.fptl} from objects of class 
\dQuote{fptl}.

\code{\link{FPTL}} to evaluate the FPTL function and create objects of class \dQuote{fptl}.
}
\examples{
## Continuing the summary.fptl(.) example:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y);
LognormalFEx <- diffproc(c("`h(t)`*x", "sigma^2*x^2", "dnorm((log(x)-(log(y)+`H(s,t)`-(sigma^2/2)*(t - s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+ `H(s,t)`-(sigma^2/2)*(t-s),sigma*sqrt(t-s))"));
z <- FPTL(dp = LognormalFEx, t0 = 1, T = 10, x0 = 1, S = 15, env = list(sigma=0.1, `h(t)` = "t/4", `H(s,t)` = "(t^2-s^2)/8")); zz <- summary(z)}
Integration.Steps(yy)
Integration.Steps(yy, from.t0 = TRUE)
Integration.Steps(yy, to.T = TRUE, n = 100, p = 0.25)

Integration.Steps(zz)
}
