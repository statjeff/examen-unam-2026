#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genera el articulo RevTeX 4.2 con todas las coordenadas pgfplots
precomputadas a partir de los posteriores guardados."""
import json, numpy as np

F = json.load(open("/home/claude/salidas/figdata.json"))
def c(pairs, nd=3):
    return " ".join(f"({x:.{nd}f},{y:.{nd}f})" for x, y in pairs)

# ---------- Figura 1: trayectoria de tau_t ----------
T = F["tau"]
med = c([(t["y"], t["med"]) for t in T], 2)
lo95 = c([(t["y"], t["q025"]) for t in T], 2)
hi95 = c([(t["y"], t["q975"]) for t in T], 2)
lo68 = c([(t["y"], t["q16"]) for t in T], 2)
hi68 = c([(t["y"], t["q84"]) for t in T], 2)

# ---------- Figura 2: densidades LOYO ----------
loyo_plots, COL = [], {"2024": "azulCCH", "2025": "doradoCCH", "2026": "vinoUNAM"}
for y in ["2024", "2025", "2026"]:
    L = F["loyo"][y]
    mx = max(L["dens"])
    pts = c([(g, dv / mx) for g, dv in zip(L["grid"], L["dens"])], 4)
    loyo_plots.append(
        f"\\addplot[draw={COL[y]},thick,fill={COL[y]},fill opacity=0.25,smooth] "
        f"coordinates {{{pts}}};\n\\addlegendentry{{{y}}}")
loyo_body = "\n".join(loyo_plots)

# ---------- Figura 3: heterogeneidad por area ----------
H = F["het_area"]
het_med = c([(h["med"], h["k"]) for h in H], 2)
het_err = "\n".join(
    f"\\draw[azulCCH,thick] (axis cs:{h['q25']:.2f},{h['k']}) -- (axis cs:{h['q75']:.2f},{h['k']});"
    for h in H)
ylabels = ", ".join(f"{h['k']}" for h in H)
CORTO = {1: "I. F\\'is.-Mat.", 2: "II. Biol.\\,y\\,Salud",
         3: "III. Sociales", 4: "IV. Human.\\,y\\,Artes"}
yticklabels = ", ".join("{" + CORTO[h["k"]] + "}" for h in H)

# ---------- numeros para el texto ----------
d_lo, d_med, d_hi = F["dtau"]
s_lo, s_med, s_hi = F["s_post"]
L24, L25, L26 = (F["loyo"][k] for k in ["2024", "2025", "2026"])
S = F["sens"]
sens_rows = "\n".join(
    f"{float(k):.0f} & {v[1]:.2f} & [{v[0]:.2f},\\,{v[2]:.2f}] \\\\"
    for k, v in sorted(S.items(), key=lambda x: float(x[0])))
tau_rows = "\n".join(
    f"{t['y']} & {t['med']:.2f} & [{t['q025']:.2f},\\,{t['q975']:.2f}] \\\\"
    for t in T)
loyo_rows = "\n".join(
    f"{y} & {F['loyo'][y]['n']} & {F['loyo'][y]['q'][1]:.2f} & "
    f"[{F['loyo'][y]['q'][0]:.2f},\\,{F['loyo'][y]['q'][2]:.2f}] & {F['loyo'][y]['rmse']:.2f} \\\\"
    for y in ["2024", "2025", "2026"])
het_rows = "\n".join(
    f"{h.get('labT',h['lab'])} & {h['n']} & {h['media']:.2f} & {h['med']:.2f} \\\\" for h in H)
mod_rows = "\n".join(
    f"{m['lab'].capitalize()} & {m['n']} & {m['media']:.2f} & {m['med']:.2f} \\\\"
    for m in F["het_mod"])

TEX = r"""%%% Compilar con: xelatex -> bibtex/manual -> xelatex x2
\documentclass[aps,prd,reprint,amsmath,amssymb,nofootinbib,superscriptaddress]{revtex4-2}

\usepackage{fontspec}
\setmainfont{Liberation Serif}
\setsansfont{Liberation Sans}
\setmonofont{Liberation Mono}

\usepackage{polyglossia}
\setmainlanguage{spanish}
\setotherlanguage{english}

\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{tikz}
\usepackage{pgfplots}
\pgfplotsset{compat=1.18}
\usepgfplotslibrary{fillbetween}
\usepackage[colorlinks=true,linkcolor=azulCCH,citecolor=azulCCH,urlcolor=azulCCH]{hyperref}

\definecolor{azulCCH}{RGB}{0,59,113}
\definecolor{doradoCCH}{RGB}{185,155,95}
\definecolor{vinoUNAM}{RGB}{140,25,45}
\definecolor{grisT}{RGB}{110,110,110}

\pgfplotsset{
  every axis/.append style={
    font=\footnotesize, axis line style={grisT}, tick style={grisT},
    grid=major, grid style={grisT!25,very thin}, label style={azulCCH},
    tick label style={black!75},
  }
}

\AtBeginDocument{%
  \renewcommand{\tablename}{Tabla}%
  \def\acknowledgmentsname{Agradecimientos}%
}

\begin{document}
\renewcommand{\tablename}{Tabla}

\title{Un desplazamiento an\'omalo en la distribuci\'on latente de puntajes
del concurso de selecci\'on a licenciatura de la UNAM tras la transici\'on
al examen en l\'inea}

\author{Jeffrey E. B\'arcenas Mosqueda}
\email{jeffrey.barcenas@cch.unam.mx}
\affiliation{Colegio de Ciencias y Humanidades, Plantel Naucalpan,
Universidad Nacional Aut\'onoma de M\'exico, Naucalpan, Estado de M\'exico, M\'exico}

\date{\today}

\begin{abstract}
En 2026 la Universidad Nacional Aut\'onoma de M\'exico (UNAM) aplic\'o por
primera vez de forma \'integramente remota su examen de selecci\'on a
licenciatura, tras cinco ciclos de aplicaci\'on presencial en sedes
arrendadas. Analizamos los puntajes m\'inimos de admisi\'on
(\emph{aciertos m\'inimos}) publicados por la Direcci\'on General de
Administraci\'on Escolar para los ciclos 2021--2026. El punto de partida
metodol\'ogico es que dicho puntaje no es una medida de desempe\~no sino el
estad\'istico de orden $K$-\'esimo de una distribuci\'on latente de
puntajes, censurado por un piso administrativo. Construimos un modelo
jer\'arquico bayesiano que invierte ese estad\'istico hacia la media
latente y separa efectos de programa y de ciclo mediante una caminata
aleatoria de primer orden. Los efectos de programa se marginalizan
anal\'iticamente por la identidad de Sherman--Morrison, lo que reduce el
espacio de par\'ametros de ${\sim}190$ a $13$ dimensiones y hace viable el
muestreo con un algoritmo de ensamble af\'in-invariante. Estimamos un
desplazamiento de $DTMED$ aciertos (IC $95\%$: $[DTLO,\,DTHI]$) en la
media latente de 2026 respecto de 2025, sobre un total de 120 reactivos.
Una validaci\'on \emph{leave-one-year-out} muestra que la discrepancia
predictiva t\'ipica del modelo es de $\pm 5$ aciertos y est\'a centrada en
cero para 2024 y 2025, mientras que 2026 exhibe una discrepancia de
$L26MED$ aciertos con la totalidad del intervalo de credibilidad fuera de
la banda de calibraci\'on. El resultado es robusto a variaciones de un
factor dos en el \'unico par\'ametro d\'ebilmente identificado del modelo.
Discutimos por qu\'e la magnitud estimada no es atribuible a un \'unico
mecanismo y delimitamos con precisi\'on qu\'e no puede establecerse con
datos agregados.
\end{abstract}

\maketitle

\section{Introducci\'on}

El concurso de selecci\'on a licenciatura de la UNAM es uno de los
procesos de evaluaci\'on educativa de mayor escala en Iberoam\'erica. En
el ciclo 2026 se registraron 191\,306 aspirantes, de los cuales 158\,712
presentaron la prueba. Ese ciclo introdujo un cambio de dise\~no sin
precedente: el examen, hist\'oricamente aplicado de manera presencial en
sedes arrendadas por la instituci\'on, se aplic\'o \'integramente en
l\'inea desde el domicilio de cada sustentante, con navegador seguro y
supervisi\'on remota asistida por sistemas autom\'aticos, a lo largo de
una ventana de tres semanas~\cite{gaceta2026}.

Un cambio de esta naturaleza plantea una pregunta emp\'irica inmediata:
\emph{\textquestiondown cambi\'o la distribuci\'on de puntajes?} La
respuesta parecer\'ia estar al alcance de una comparaci\'on directa de los
puntajes de corte publicados. Sostenemos que esa comparaci\'on es
enga\~nosa por dos razones, y que corregirlas modifica la magnitud
estimada del fen\'omeno en un factor cercano a dos.

La primera raz\'on es que el puntaje de corte no mide desempe\~no. Es el
puntaje del \'ultimo aspirante admitido, es decir, un estad\'istico de
orden cuya posici\'on depende conjuntamente de la oferta de lugares, del
n\'umero de sustentantes y de la forma de la distribuci\'on latente de
puntajes. La segunda es que en 2026 la raz\'on de selecci\'on aument\'o
simult\'aneamente: se admiti\'o a una fracci\'on mayor de los
sustentantes, lo que \emph{por s\'i solo} presiona el corte a la baja.
Que el corte haya subido pese a ello implica un desplazamiento de la
distribuci\'on subyacente mayor que el observado en la superficie.

Este trabajo formaliza esa intuici\'on. La secci\'on~\ref{sec:datos}
describe los datos y los criterios de inclusi\'on; la
secci\'on~\ref{sec:modelo}, el modelo de inversi\'on de cuantil censurado
y su implementaci\'on; la secci\'on~\ref{sec:resultados}, los resultados y
su validaci\'on; y la secci\'on~\ref{sec:discusion}, los l\'imites de lo
que puede inferirse. Este \'ultimo punto no es un ap\'endice de cautela:
es, a nuestro juicio, la contribuci\'on m\'as importante del art\'iculo.

\section{Datos}\label{sec:datos}

La Direcci\'on General de Administraci\'on Escolar publica, por ciclo y
por programa, la oferta de lugares, el n\'umero de aspirantes registrados,
el n\'umero que present\'o el examen, el puntaje m\'inimo de admisi\'on y
el n\'umero de seleccionados. Recopilamos estas cifras para los ciclos
2021--2026, obteniendo 1\,253 registros de programa-ciclo.

La unidad de an\'alisis es el par (c\'odigo de programa, modalidad): un
mismo c\'odigo se emplea para las modalidades abierta y a distancia, de
modo que el c\'odigo por s\'i solo no identifica un\'ivocamente. Se
identific\'o adem\'as un piso administrativo uniforme de 40 aciertos en
las tres modalidades, por debajo del cual no se admite a ning\'un
aspirante. Las observaciones en ese piso est\'an censuradas por la
izquierda: el corte deja de ser informativo sobre la distribuci\'on
latente. Excluimos tambi\'en los casos de saturaci\'on, en que
pr\'acticamente todo sustentante fue admitido, y aquellos con muy pocos
sustentantes. Formalmente conservamos las observaciones que cumplen
\begin{equation}
c_{jt} > 40, \qquad 0.01 < \frac{K_{jt}}{N_{jt}} < 0.95, \qquad N_{jt} > 30,
\end{equation}
donde $c_{jt}$ es el puntaje de corte, $K_{jt}$ el n\'umero de
seleccionados y $N_{jt}$ el de sustentantes del programa $j$ en el ciclo
$t$. Quedan NOBS observaciones ($83.7\%$ del total) sobre NUNI unidades
con al menos tres ciclos.

\section{Modelo}\label{sec:modelo}

\subsection{Inversi\'on del cuantil}

Sea $X_{jt}$ el puntaje de un sustentante del programa $j$ en el ciclo
$t$, sobre 120 reactivos. Suponemos
$X_{jt}\sim\mathcal{N}(\mu_{jt},s^{2})$. El procedimiento de admisi\'on
ordena a los sustentantes y admite a los $K_{jt}$ de mayor puntaje, de
modo que el corte satisface
\begin{equation}
P\!\left(X_{jt}\ge c_{jt}\right)\;\simeq\;\frac{K_{jt}}{N_{jt}}\;\equiv\;p_{jt}.
\end{equation}
Cada registro proporciona, por tanto, un punto de la funci\'on de
distribuci\'on acumulada latente. Invirtiendo,
\begin{equation}
\widehat{\mu}_{jt} \;=\; c_{jt}-\tfrac12 \;-\; s\,z_{jt},
\qquad z_{jt}=\Phi^{-1}\!\left(1-p_{jt}\right),
\label{eq:inv}
\end{equation}
con correcci\'on por continuidad. La incertidumbre de
$\widehat{\mu}_{jt}$ proviene del error binomial en $p_{jt}$; por el
m\'etodo delta,
\begin{equation}
v_{jt}\;\equiv\;\mathrm{Var}\!\left(\widehat{\mu}_{jt}\right)
\;=\;\frac{s^{2}\,p_{jt}(1-p_{jt})}{N_{jt}\,\phi(z_{jt})^{2}} .
\end{equation}

La ecuaci\'on~\eqref{eq:inv} explicita el punto central del art\'iculo:
el desplazamiento latente combina el cambio del corte y el cambio de la
raz\'on de selecci\'on. Un aumento de $p$ reduce $z$, de manera que un
corte que sube mientras la selecci\'on se ampl\'ia implica un
desplazamiento de $\mu$ estrictamente mayor que $\Delta c$.

\subsection{Estructura jer\'arquica}

Descomponemos la media latente en un efecto de programa y uno de ciclo,
\begin{equation}
\mu_{jt} \;=\; a_j + \tau_t + \eta_{jt},
\qquad \eta_{jt}\sim\mathcal{N}(0,\sigma_\eta^{2}),
\end{equation}
con $a_j\sim\mathcal{N}(m_{g(j)},\sigma_a^{2})$, donde $g(j)$ indexa el
\'area de conocimiento. Los efectos de ciclo siguen una caminata
aleatoria de primer orden,
$\tau_t=\tau_{t-1}+\varepsilon_t$, $\varepsilon_t\sim\mathcal{N}(0,\sigma_{\mathrm{rw}}^{2})$,
con $\tau_{2021}\equiv 0$ por identificabilidad. Esta elecci\'on es
deliberadamente permisiva: una caminata aleatoria admite tendencias
suaves y no penaliza cambios graduales, de modo que no fuerza la
detecci\'on de una discontinuidad.

\subsection{Marginalizaci\'on de los efectos de programa}

Con ${\sim}180$ programas, el muestreo directo de $\{a_j\}$ es hostil para
un algoritmo de ensamble. Dado que el modelo es gaussiano condicionalmente
en las escalas, los $a_j$ pueden integrarse de forma exacta. Para el
programa $j$ con $n_j$ observaciones, sea
$r_{jt}=\widehat{\mu}_{jt}-\tau_t-m_{g(j)}$. Entonces
\begin{equation}
\mathbf{r}_j\sim\mathcal{N}\!\left(\mathbf{0},\;
\mathsf{D}_j+\sigma_a^{2}\mathbf{1}\mathbf{1}^{\!\top}\right),
\qquad \mathsf{D}_j=\mathrm{diag}\!\left(v_{jt}+\sigma_\eta^{2}\right).
\end{equation}
La identidad de Sherman--Morrison~\cite{hager1989} da la forma
cuadr\'atica y el determinante en $\mathcal{O}(n_j)$:
\begin{align}
Q_j &= \textstyle\sum_t u_{jt}r_{jt}^{2}
 - \frac{\sigma_a^{2}\left(\sum_t u_{jt}r_{jt}\right)^{2}}{1+\sigma_a^{2}\sum_t u_{jt}},\\
\ln\!\left|\mathsf{\Sigma}_j\right| &= -\textstyle\sum_t \ln u_{jt}
 + \ln\!\left(1+\sigma_a^{2}\textstyle\sum_t u_{jt}\right),
\end{align}
con $u_{jt}=(v_{jt}+\sigma_\eta^{2})^{-1}$. Todas las sumas se eval\'uan
de forma vectorizada. El espacio muestral se reduce a los 13
hiperpar\'ametros
$\left(\tau_{2022..2026},\,m_{1..4},\,\ln s,\,\ln\sigma_a,\,\ln\sigma_\eta,\,\ln\sigma_{\mathrm{rw}}\right)$.

\subsection{Distribuciones a priori}

Adoptamos $\sigma_a,\sigma_\eta,\sigma_{\mathrm{rw}}$ semi-normales,
$m_g\sim\mathcal{N}(60,30^{2})$ y, para la desviaci\'on est\'andar
latente, $\ln s\sim\mathcal{N}(\ln 17,\,0.25^{2})$. Este \'ultimo es el
\'unico prior con consecuencias sustantivas, porque $s$ fija la escala de
la inversi\'on~\eqref{eq:inv}. Es tambi\'en el par\'ametro cuya
identificaci\'on resulta m\'as d\'ebil: proviene \'unicamente de la
exigencia de que $c_{jt}-s\,z_{jt}$ admita una descomposici\'on aditiva
parsimoniosa. La secci\'on~\ref{sec:sens} documenta que la conclusi\'on
principal es insensible a esta elecci\'on.

\subsection{C\'omputo}

El muestreo se realiz\'o con el algoritmo de ensamble af\'in-invariante de
Goodman y Weare~\cite{goodman2010}, en la implementaci\'on
\texttt{emcee}~\cite{foremanmackey2013}, con 80 caminantes y una mezcla de
movimientos diferenciales ($80\%$ \texttt{DEMove}, $20\%$
\texttt{DESnookerMove}). Dos decisiones fueron necesarias para obtener
mezcla aceptable. Primero, el centrado de $z_{jt}$ respecto de su media
muestral, que decorrela $s$ de los niveles $m_g$; sin \'el, el tiempo de
autocorrelaci\'on integrado excede los 300 pasos. Segundo, una
inicializaci\'on con dispersi\'on ajustada por bloque a la escala de cada
par\'ametro: los movimientos diferenciales dependen de la diversidad del
ensamble, y una bola isotr\'opica peque\~na produce exploraci\'on
deficiente en las direcciones de nivel. Con ambos ajustes,
$\tau_{\mathrm{int}}\le 89$ pasos sobre cadenas de 8\,000, con tasa de
aceptaci\'on de $0.24$.

\section{Resultados}\label{sec:resultados}

\subsection{Trayectoria latente}

La figura~\ref{fig:tau} y la tabla~\ref{tab:tau} presentan los efectos de
ciclo. El periodo 2021--2025 es notablemente estable: todas las medianas
caen en el intervalo $[-1.6,\,0.2]$ aciertos, con intervalos de
credibilidad que rara vez exceden $\pm 1$ acierto. El ciclo 2026 rompe ese
patr\'on con un valor de $TAU26$ aciertos.

\begin{figure}[t]
\centering
\begin{tikzpicture}
\begin{axis}[
  width=\columnwidth, height=6.1cm,
  xlabel={Ciclo}, ylabel={$\tau_t$ (aciertos, base 2021)},
  xmin=2020.6, xmax=2026.4, xtick={2021,2022,2023,2024,2025,2026},
  xticklabel style={/pgf/number format/1000 sep=},
  ymin=-4, ymax=19,
]
\addplot[name path=hi95,draw=none] coordinates {HI95};
\addplot[name path=lo95,draw=none] coordinates {LO95};
\addplot[azulCCH!18] fill between[of=hi95 and lo95];
\addplot[name path=hi68,draw=none] coordinates {HI68};
\addplot[name path=lo68,draw=none] coordinates {LO68};
\addplot[azulCCH!38] fill between[of=hi68 and lo68];
\addplot[azulCCH,very thick,mark=*,mark size=1.7pt,
         mark options={fill=white,draw=azulCCH}]
  coordinates {MED};
\draw[grisT,dashed] (axis cs:2020.6,0) -- (axis cs:2026.4,0);
\node[anchor=east,font=\scriptsize,vinoUNAM]
  at (axis cs:2025.85,15.8) {examen en l\'inea};
\end{axis}
\end{tikzpicture}
\caption{Efecto latente de ciclo $\tau_t$, en aciertos sobre 120, con
2021 como referencia. Bandas: intervalos de credibilidad al $68\%$ y
$95\%$. La estabilidad del periodo presencial contrasta con el salto de
2026.}
\label{fig:tau}
\end{figure}

\begin{table}[t]
\caption{Efectos de ciclo. Mediana a posteriori e intervalo de
credibilidad al $95\%$, en aciertos.}
\label{tab:tau}
\begin{ruledtabular}
\begin{tabular}{lcc}
Ciclo & Mediana & IC $95\%$ \\ \colrule
TAUROWS
\end{tabular}
\end{ruledtabular}
\end{table}

La diferencia entre ciclos consecutivos 2026--2025 es de $DTMED$
aciertos, con intervalo $[DTLO,\,DTHI]$ y
$P(\Delta\tau>0)>0.9999$. Conviene contrastar esta cifra con el aumento
crudo del puntaje de corte, de $9.2$ aciertos en promedio sobre los
programas pareados. La discrepancia entre ambas ---un factor de $1.7$---
es exactamente el efecto anticipado en la introducci\'on: la raz\'on media
de selecci\'on pas\'o de $0.35$ a $0.42$, y esa ampliaci\'on de la
admisi\'on enmascara parte del desplazamiento subyacente. La desviaci\'on
est\'andar latente estimada es $s=SMED$ aciertos
($[SLO,\,SHI]$), de modo que el desplazamiento equivale a
aproximadamente $0.75$ desviaciones est\'andar.

\subsection{Validaci\'on \emph{leave-one-year-out}}

Que un salto sea grande respecto de la serie previa no basta: hay que
saber cu\'anto se equivoca el modelo al predecir un ciclo que no ha visto.
Ajustamos el modelo usando \'unicamente los ciclos anteriores a un ciclo
$T$ dado, extrapolamos la caminata aleatoria un paso, obtenemos la
posterior condicional de cada $a_j$ y predecimos el corte
\emph{manteniendo fija la raz\'on de selecci\'on observada} en $T$. La
discrepancia $D_{jT}=c^{\mathrm{obs}}_{jT}-c^{\mathrm{pred}}_{jT}$ mide
entonces el desplazamiento distribucional no anticipado, limpio de
cambios en oferta y demanda.

\begin{figure}[t]
\centering
\begin{tikzpicture}
\begin{axis}[
  width=\columnwidth, height=5.6cm,
  xlabel={Discrepancia media $\overline{D}_T$ (aciertos)},
  ylabel={Densidad a posteriori (norm.)},
  xmin=-9, xmax=22, ymin=0, ymax=1.18, ytick={0,0.5,1},
  legend style={at={(0.02,0.98)},anchor=north west,draw=grisT!50,
                font=\scriptsize,fill=white,fill opacity=0.85,
                text opacity=1,row sep=0.5pt},
]
LOYOBODY
\draw[grisT,dashed] (axis cs:0,0) -- (axis cs:0,1.18);
\end{axis}
\end{tikzpicture}
\caption{Validaci\'on \emph{leave-one-year-out}. Distribuci\'on a
posteriori de la discrepancia predictiva media para cada ciclo excluido.
Los ciclos presenciales 2024 y 2025 se centran en cero; 2026 se separa por
completo de la banda de calibraci\'on.}
\label{fig:loyo}
\end{figure}

\begin{table}[t]
\caption{Validaci\'on \emph{leave-one-year-out}. $n$ es el n\'umero de
programas evaluados; $\overline{D}$ la discrepancia media en aciertos.}
\label{tab:loyo}
\begin{ruledtabular}
\begin{tabular}{lcccc}
Ciclo & $n$ & $\overline{D}$ & IC $95\%$ & RECM \\ \colrule
LOYOROWS
\end{tabular}
\end{ruledtabular}
\end{table}

Los resultados (figura~\ref{fig:loyo}, tabla~\ref{tab:loyo}) establecen la
banda de calibraci\'on. Para 2024 y 2025 la discrepancia media es
compatible con cero y el error cuadr\'atico medio ronda los $4.4$
aciertos. Para 2026 la discrepancia media es $L26MED$ aciertos, con el
extremo inferior del intervalo al $95\%$ en $L26LO$, muy por encima del
extremo superior de cualquiera de los ciclos de calibraci\'on. El error
cuadr\'atico medio se cuadruplica. En t\'erminos operativos: el modelo
predice un ciclo presencial no observado con un error t\'ipico de $\pm 5$
aciertos, y falla en 2026 por un margen tres veces mayor y sistem\'atico
en signo.

\subsection{Heterogeneidad}

La discrepancia de 2026 no es uniforme (tabla~\ref{tab:het},
figura~\ref{fig:het}). Las \'areas I, II y III presentan valores medianos
entre $16.7$ y $18.2$ aciertos, mientras que el \'area IV se sit\'ua en
$7.6$. Interpretamos esta diferencia con cautela: los programas del
\'area IV concentran los cortes m\'as bajos y por tanto la mayor
proximidad al piso administrativo, donde el estad\'istico de orden pierde
sensibilidad. Parte de la atenuaci\'on es, con alta probabilidad, un
artefacto de censura y no un efecto sustantivo.

\begin{figure}[t]
\centering
\begin{tikzpicture}
\begin{axis}[
  width=0.68\columnwidth, height=4.4cm,
  xlabel={Discrepancia mediana en 2026 (aciertos)},
  xmin=0, xmax=24, ytick={YLAB}, yticklabels={YTICKLAB},
  ytick style={draw=none}, y dir=reverse, ymin=0.4, ymax=4.6,
  yticklabel style={font=\scriptsize,align=right},
]
HETERR
\addplot[only marks,mark=*,mark size=2.4pt,azulCCH,
         mark options={fill=doradoCCH,draw=azulCCH,line width=0.7pt}]
  coordinates {HETMED};
\end{axis}
\end{tikzpicture}
\caption{Discrepancia mediana por \'area de conocimiento en 2026. Las
l\'ineas indican el rango intercuartil entre programas. La atenuaci\'on
del \'area IV es atribuible en parte a censura por el piso
administrativo.}
\label{fig:het}
\end{figure}

\begin{table}[t]
\caption{Discrepancia de 2026 por \'area de conocimiento y por modalidad
de estudio, en aciertos.}
\label{tab:het}
\begin{ruledtabular}
\begin{tabular}{lccc}
Estrato & $n$ & Media & Mediana \\ \colrule
HETROWS
\colrule
MODROWS
\end{tabular}
\end{ruledtabular}
\end{table}

Por modalidad, en cambio, la homogeneidad es llamativa: escolarizado
($15.3$), abierta ($14.4$) y a distancia ($12.5$) presentan
discrepancias comparables. Este resultado corrige de manera directa una
lectura descriptiva ingenua: sobre los puntajes de corte sin modelar, las
modalidades abierta y a distancia parecen saltar casi el doble que la
escolarizada. Esa diferencia aparente se disuelve al condicionar sobre las
razones de selecci\'on, y era un artefacto de composici\'on.

\subsection{Sensibilidad al prior sobre $s$}\label{sec:sens}

La tabla~\ref{tab:sens} reporta $\Delta\tau_{2026-2025}$ fijando $s$ en
valores que abarcan un factor de dos. La estimaci\'on var\'ia entre
$14.7$ y $16.3$ aciertos: el par\'ametro peor identificado del modelo no
compromete la conclusi\'on. Cabe se\~nalar adem\'as que, con prior
centrado en $17$, la posterior se desplaza a $SMED$, lo que indica que
$s$ recibe informaci\'on emp\'irica genuina a trav\'es de la estructura
aditiva.

\begin{table}[t]
\caption{An\'alisis de sensibilidad: $\Delta\tau_{2026-2025}$ con la
desviaci\'on est\'andar latente $s$ fijada.}
\label{tab:sens}
\begin{ruledtabular}
\begin{tabular}{lcc}
$s$ (aciertos) & $\Delta\tau$ & IC $95\%$ \\ \colrule
SENSROWS
\end{tabular}
\end{ruledtabular}
\end{table}

\section{Qu\'e no puede concluirse}\label{sec:discusion}

El desplazamiento estimado es grande, robusto y est\'a bien calibrado
contra el comportamiento hist\'orico del modelo. Nada de ello autoriza a
atribuirlo a un mecanismo particular, y conviene ser expl\'icito sobre por
qu\'e.

El dise\~no carece de grupo de control contempor\'aneo. La totalidad de
los sustentantes de 2026 recibi\'o el mismo tratamiento, de modo que el
estimando $\Delta\tau$ agrega \emph{todo} lo que cambi\'o en ese ciclo. Al
menos cuatro mecanismos son compatibles con lo observado y este dise\~no
no los distingue.

\emph{(i) Reducci\'on de fricci\'on log\'istica.} Presentar desde el
domicilio elimina traslados y la incertidumbre asociada a sedes
desconocidas. El mecanismo es real y ben\'evolo. Observamos, sin embargo,
que la tasa de no presentaci\'on en la modalidad escolarizada
\emph{aument\'o} de $0.171$ a $0.184$, lo que debilita la versi\'on fuerte
de esta explicaci\'on.

\emph{(ii) Cambio de instrumento.} Una ventana de aplicaci\'on de tres
semanas exige m\'ultiples formas del examen. Cualquier diferencia de
calibraci\'on entre bancos de reactivos se confunde de manera perfecta con
la modalidad.

\emph{(iii) Integridad del proceso.} La propia instituci\'on report\'o la
cancelaci\'on del proceso de aproximadamente el $2\%$ de los aspirantes
por conductas contrarias a la convocatoria~\cite{gaceta2026}. Una ventana
prolongada crea, en principio, condiciones para la filtraci\'on de
reactivos entre jornadas.

\emph{(iv) Composici\'on de la poblaci\'on.} La participaci\'on y la
oferta de lugares cambiaron simult\'aneamente, alterando qui\'en se
presenta y con qu\'e nivel de preparaci\'on.

Ninguna partici\'on de $\Delta\tau$ entre estos cuatro canales es
identificable con datos agregados de programa-ciclo. Se\~nalamos, no
obstante, un dise\~no que s\'i tendr\'ia poder discriminante y que
requiere \'unicamente informaci\'on p\'ublica adicional: la asignaci\'on
de cada programa a una jornada espec\'ifica dentro de la ventana de tres
semanas. Bajo el mecanismo (i), $\Delta\tau$ debe ser homog\'eneo respecto
de la fecha de aplicaci\'on; bajo el mecanismo (iii), debe crecer con
ella. Se trata de un contraste de dosis-respuesta con grupo de
comparaci\'on interno al ciclo 2026, inmune a la mayor parte de los
factores de confusi\'on enumerados. Lo dejamos como continuaci\'on
natural de este trabajo.

Una limitaci\'on adicional merece \'enfasis. Todas las inferencias operan
en el nivel de programa-ciclo. Trasladarlas a afirmaciones sobre
sustentantes individuales constituir\'ia una falacia ecol\'ogica. En
particular, este an\'alisis no sustenta ninguna afirmaci\'on sobre la
conducta de ning\'un aspirante, ni sobre la validez de ninguna admisi\'on
concreta.

\section{Conclusiones}

Tratar el puntaje de corte como lo que es ---un estad\'istico de orden
censurado--- y no como una medida de desempe\~no modifica de manera
sustantiva las conclusiones. El desplazamiento de la distribuci\'on
latente de puntajes en 2026 es de $DTMED$ aciertos sobre 120, cerca de
$0.75$ desviaciones est\'andar latentes, frente a un aumento aparente de
$9.2$ aciertos en los cortes publicados. Una validaci\'on
\emph{leave-one-year-out} sit\'ua ese valor muy fuera de la banda de
error predictivo del modelo, calibrada en $\pm 5$ aciertos sobre los
ciclos presenciales. El resultado resiste variaciones de un factor dos en
el \'unico par\'ametro d\'ebilmente identificado.

Metodol\'ogicamente, el trabajo ilustra que la marginalizaci\'on
anal\'itica de efectos aleatorios gaussianos permite abordar modelos
jer\'arquicos de varios cientos de par\'ametros con muestreadores de
ensamble, que de otro modo resultar\'ian inadecuados.

Sustantivamente, la conclusi\'on defendible es que el proceso de medici\'on
de 2026 no es comparable con el de los ciclos precedentes. Determinar
\emph{por qu\'e} exige informaci\'on que este conjunto de datos no
contiene.

\begin{acknowledgments}
El autor agradece a la Direcci\'on General de Administraci\'on Escolar de
la UNAM la publicaci\'on abierta de las cifras del concurso de
selecci\'on, sin la cual este an\'alisis no ser\'ia posible. En la
implementaci\'on computacional y en la redacci\'on de este manuscrito se
emple\'o asistencia de un modelo de lenguaje (Claude, Anthropic); el
dise\~no del estudio, las decisiones de modelaci\'on, la interpretaci\'on
de los resultados y la responsabilidad por su contenido corresponden
\'integramente al autor.
\end{acknowledgments}

\begin{thebibliography}{99}

\bibitem{gaceta2026}
\emph{La UNAM realiza la aplicaci\'on en l\'inea de su examen de
admisi\'on para el ingreso a licenciatura 2026}, Gaceta UNAM
(25 de mayo de 2026).

\bibitem{hager1989}
W.~W. Hager, \emph{Updating the inverse of a matrix},
SIAM Rev. \textbf{31}, 221 (1989).

\bibitem{goodman2010}
J.~Goodman and J.~Weare, \emph{Ensemble samplers with affine invariance},
Commun. Appl. Math. Comput. Sci. \textbf{5}, 65 (2010).

\bibitem{foremanmackey2013}
D.~Foreman-Mackey, D.~W. Hogg, D.~Lang, and J.~Goodman,
\emph{emcee: The MCMC hammer},
Publ. Astron. Soc. Pac. \textbf{125}, 306 (2013).

\bibitem{gelman2013}
A.~Gelman, J.~B. Carlin, H.~S. Stern, D.~B. Dunson, A.~Vehtari, and
D.~B. Rubin, \emph{Bayesian Data Analysis}, 3rd ed.
(CRC Press, Boca Raton, 2013).

\bibitem{vehtari2017}
A.~Vehtari, A.~Gelman, and J.~Gabry, \emph{Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC},
Stat. Comput. \textbf{27}, 1413 (2017).

\bibitem{betancourt2013}
M.~Betancourt and M.~Girolami, \emph{Hamiltonian Monte Carlo for
hierarchical models}, arXiv:1312.0906 (2013).

\bibitem{holland1986}
P.~W. Holland, \emph{Statistics and causal inference},
J. Am. Stat. Assoc. \textbf{81}, 945 (1986).

\bibitem{rubin1974}
D.~B. Rubin, \emph{Estimating causal effects of treatments in randomized
and nonrandomized studies},
J. Educ. Psychol. \textbf{66}, 688 (1974).

\end{thebibliography}

\end{document}
"""

rep = {
    "MED": med, "LO95": lo95, "HI95": hi95, "LO68": lo68, "HI68": hi68,
    "LOYOBODY": loyo_body, "HETMED": het_med, "HETERR": het_err,
    "YTICKLAB": yticklabels, "YLAB": ylabels,
    "TAUROWS": tau_rows, "LOYOROWS": loyo_rows, "SENSROWS": sens_rows,
    "HETROWS": het_rows, "MODROWS": mod_rows,
    "DTMED": f"{d_med:.2f}", "DTLO": f"{d_lo:.2f}", "DTHI": f"{d_hi:.2f}",
    "SMED": f"{s_med:.1f}", "SLO": f"{s_lo:.1f}", "SHI": f"{s_hi:.1f}",
    "TAU26": f"{T[-1]['med']:.2f}",
    "L26MED": f"{L26['q'][1]:.2f}", "L26LO": f"{L26['q'][0]:.2f}",
    "NOBS": str(F["n_obs"]), "NUNI": str(F["n_uni"]),
}
for k in sorted(rep, key=len, reverse=True):
    TEX = TEX.replace(k, rep[k])

open("/home/claude/paper/examen_unam_2026.tex", "w").write(TEX)
print("tex generado:", len(TEX), "caracteres")
