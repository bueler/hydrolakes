\section{Alternative: water pressure is a state variable} \label{app:alternative}  Our model with state variables $(W,Y)$ above, concretely equations \eqref{eq:Weqnforsemi} and \eqref{eq:Yeqnforsemi}, has disadvantages.  We now describe an alternative formulation which has state variables $(W,P)$.  Here is a sketch of the relative advantages of a $(W,P)$ model:\begin{itemize}
\item the pressure $P$ is more directly-related to ice dynamics than is $Y$, because of equation \eqref{eq:mohr-coulomb},
\item the variable $Y$ is not more observable than $P$, and it may be the opposite because $P$ can be observed in bore holes,
\item because the new model has a time-derivative expression for $P$, we may be able to avoid the technique in the next appendix if stiffness arises as an issue,
\item the model in \cite{Schoofmeltsupply} is in terms of effective pressure, not water amount, and thereby the equation by which $P$ is evolving here can be more directly compared to the theory in \cite{Schoofmeltsupply},
\item both $W$ and $P$ are already in use in PISM as state variables (or nearly so).
\end{itemize}

We start by recalling the relation \eqref{eq:pressure} between $(W,Y)$ and $P$:
\begin{equation}\label{eq:pressureagain}
   P = p_i \left(\frac{W}{Y+Y_{min}}\right)^\sigma.
\end{equation}
The basic idea for converting from a $(W,Y)$ model to a $(W,P)$ model is to compute logarithms of \eqref{eq:pressureagain},
  $$\ln P = \ln p_i + \sigma (\ln W - \ln (Y+Y_{min})),$$
and then find time derivatives of both sides.  The term for the rate of change of the logarithm of overburden pressure is now assumed to be negligible, $\partial (\ln p_i)/\partial t \approx 0$, because changes to $p_i$ are on ice-dynamical timescales which are slower than the timescales from changes to water amount $W$, or even the capacity variable $Y$.  In fact, however, keeping this term requires no essential changes to our model.

Thus, recalling  $d(\ln f)/dt = (1/f) df/dt$,
\begin{equation}\label{eq:dtrelation}
   \frac{\partial P}{\partial t} = \sigma P \left(\frac{\partial(\ln W)}{\partial t} - \frac{\partial(\ln (Y+Y_{min}))}{\partial t}\right).
\end{equation}
From \eqref{eq:Weqnforsemi} and \eqref{eq:Yeqnforsemi}, respectively, we get
\begin{align}
&\frac{\partial (\ln W)}{\partial t} = \frac{1}{W} \frac{\partial W}{\partial t} = \frac{c_1}{W} \nabla \cdot \big(W \nabla P\big) + c_2 \left|\nabla P\right|^2 + \frac{S}{\rho_w W},  \label{eq:lnWeqnalt} \\
&\frac{\partial (\ln (Y+Y_{min}))}{\partial t} = \frac{1}{Y+Y_{min}} \frac{\partial Y}{\partial t}  \label{eq:lnYeqnalt} \\
   &\qquad\qquad = c_3 \left(\frac{W}{Y+Y_{min}}\right) \left|\nabla P\right|^2 - \Creep A \max\left\{0,\big(p_i - P\big)^n\right\} \frac{Y}{Y+Y_{min}}. \notag
\end{align}

We seek to write these equations without ``$Y$'' on their right sides, because we want an equation for $\partial P/\partial t$ in terms of $W$ and $P$ but not $Y$.  First note that from \eqref{eq:pressureagain}, 
	$$\frac{W}{Y+Y_{min}} = \left(\frac{P}{p_i}\right)^{1/\sigma}.$$
Second, 
	$$\frac{Y}{Y+Y_{min}} = 1 - \frac{Y_{min}}{Y+Y_{min}} = 1 - \left(\frac{P}{p_i}\right)^{1/\sigma} \frac{Y_{min}}{W},$$
also from \eqref{eq:pressureagain}.  Thus we can write out an equation for $\partial P/\partial t$ by combining \eqref{eq:dtrelation}, \eqref{eq:lnWeqnalt}, and \eqref{eq:lnYeqnalt} and these replacements.

We now have a model with state variables $(W,P)$.  But we add one more regularization, replacing the division by $W$ in equation \eqref{eq:lnWeqnalt}, which corresponds to an unbounded term as $W\to 0$, with division by $W+W_{min}$ where $W_{min}>0$ is constant and much smaller than one meter.  Now we write the $(W,P)$ form:
\begin{align}
\frac{\partial W}{\partial t} &= c_1 \nabla \cdot \big(W \nabla P\big) + c_2 W \left|\nabla P\right|^2 + \frac{S}{\rho_w},  \label{eq:Weqnalt} \\
\frac{\partial P}{\partial t} &= \sigma \frac{P}{W+W_{min}} \bigg[c_1 \nabla \cdot \big(W \nabla P\big) + \bigg(c_2 - c_3 \left(\frac{P}{p_i}\right)^{1/\sigma}\bigg) W \left|\nabla P\right|^2 \label{eq:Peqnalt} \\
     &\qquad\qquad + \Creep A \max\left\{0,\big(p_i - P\big)\right\}^n \left(W - \left(\frac{P}{p_i}\right)^{1/\sigma} Y_{min}\right) + \frac{S}{\rho_w} \bigg]. \notag 
\end{align}
The same model can be written as a system of two equations with a nondegenerate (if $W_{min}>0$) triangular form on the left:
$$\begin{bmatrix}
     1 & 0 \\ -\sigma P\phantom{\Big|} & (W+W_{min})
  \end{bmatrix}
  \begin{bmatrix}
     \frac{\partial W}{\partial t} \\ \frac{\partial P}{\partial t}\phantom{\Big|}
  \end{bmatrix} = 
  \begin{bmatrix}
     c_1 \nabla \cdot \big(W \nabla P\big) + c_2 W \left|\nabla P\right|^2 + \rho_w^{-1} S \\
     \sigma p_i^{-1/\sigma} P\, \mathcal{R}
  \end{bmatrix}$$
where $\mathcal{R}=\mathcal{R}(W,P,\nabla P)$ is this expression
	$$\mathcal{R} = - c_3 P^{1/\sigma} W \left|\nabla P\right|^2  + \Creep A \max\left\{0,\big(p_i - P\big)\right\}^n \left(p_i^{1/\sigma} W - P^{1/\sigma} Y_{min}\right),$$
which only involves the opening and closure processes.

FIXME: we can numerically solve equations \eqref{eq:Weqnalt} and \eqref{eq:Peqnalt}, implicitly and with Newton, but non-physical instabilities arise; this seems to occur because $\partial W/\partial t$ depends on spatial second derivatives of $P$, but $P$ does not evolve in a strongly smoothing way; indeed it essentially evolves as a highly nonlinear first-order equation

We must reconsider the verification we accomplished in section \ref{sec:verif}.  If we set opening and closing processes to zero ($c_2=0$, $c_3=0$, $C_{reep}=0$) and we have no water source ($S=0$) and we remove the last regularization ($W_{min}=0$) then equations \eqref{eq:Weqnalt} and \eqref{eq:Peqnalt} reduce to simpler forms:
\begin{equation}\label{eq:WPforverif}
\frac{\partial W}{\partial t} = c_1 \nabla \cdot \big(W \nabla P\big), \qquad \frac{\partial P}{\partial t} = \sigma P \,\frac{c_1}{W} \nabla \cdot \big(W \nabla P\big).
\end{equation}
Multiplying the first by $\sigma P$ and the second by $W$ and subtracting gives
  $$\sigma P \frac{\partial W}{\partial t} - W \frac{\partial P}{\partial t} = 0$$
or, equivalently, $\partial/\partial t\left(W^{-\sigma} P\right) = 0$.  That is, in the case where there is no opening or closing or sources then $W^{-\sigma} P$ is constant in time.  But of course we want $P = p_i (W/(Y+Y_{min}))^{\sigma}$, so we conclude $Y+Y_{min}$ is independent of time.  Thus the capacity thickness $Y=Y(x,y)$ is independent of time, while on the other hand the lack of opening and closing processes would imply from equation \eqref{eq:capacityevolution} that $Y$ was independent of time also.  Replacing $P$ in the first of equations \eqref{eq:WPforverif} by $P=p_i(W/Y_0)^\sigma$, which chooses a specific capacity thickness, and with $p_i = \rho_i g H_0$ which chooses a specific overburden distribution, we recover the vertically-integrated porous medium equation \eqref{eq:vertporous} of section \ref{sec:verif}, namely
	$$\frac{\partial W}{\partial t} = \gamma \nabla \cdot \big(W \nabla (W^{\sigma})\big).$$
This is the equation which has a Barenblatt exact solution, so that exact solution can be used to verify a numerical solver for coupled equations \eqref{eq:Weqnalt} and \eqref{eq:Peqnalt}.

FIXME: address steady state: if we look at steady state for equations \eqref{eq:Weqnalt} and \eqref{eq:Peqnalt} and if we consider the natural case where $p_i^{1/\sigma} W \gg P^{1/\sigma} Y_{min}$ then we get a steady state equation without ``$W$'', which should determine the distribution of $P$:
	$$|\nabla P|^2 - \frac{C_{reep} A}{c_3} \frac{p_i^{1/\sigma} (p_i-P)^n}{P^{1/\sigma}} = 0$$
this is a Hamilton-Jacobi equation $|\nabla P| = f(P)$ where $f$ is badly-behaved at $P=0$; presumably this equation applies over the set $\{W>0\}$, and if we know $P$ at the boundary of this set then apparently we can determine $P$ without knowing anything else about $W$

FIXME: we could address verification by exact solution of the steady equation for $P$, if we had that
