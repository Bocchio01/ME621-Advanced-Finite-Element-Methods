% assignment.cls
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{assignment-cls}[2024/01/19 v1.0 ME621 Assignments class]

\LoadClass{article}

% Page layout
\RequirePackage[a4paper, margin=20mm]{geometry}

% Colors and graphics
\RequirePackage[dvipsnames,table]{xcolor}
\RequirePackage{graphicx}
\RequirePackage{subcaption}

% Language settings
\RequirePackage[english]{babel}

% Hyperlinks
\RequirePackage[hidelinks]{hyperref}

% Mathematical typesetting
\RequirePackage{amsmath}

% Code listings
\RequirePackage{listings}
\RequirePackage{xcolor}

\definecolor{codegreen}{rgb}{0,0.6,0}
\definecolor{codegray}{rgb}{0.5,0.5,0.5}
\definecolor{codepurple}{rgb}{0.58,0,0.82}
\definecolor{backcolour}{rgb}{0.95,0.95,0.92}

\lstset{
    basicstyle          =   \sffamily,
    keywordstyle        =   \bfseries,
    commentstyle        =   \rmfamily\itshape,
    stringstyle         =   \ttfamily,
    flexiblecolumns,
    numbers             =   left,
    showspaces          =   false,
    showstringspaces    =   false,
    captionpos          =   t,
    frame               =   single,
    xleftmargin=15pt,
    xrightmargin=5pt,
}

\lstdefinestyle{mathematica}{
    backgroundcolor=\color{white},
    commentstyle=\color[HTML]{338AAF}\ttfamily,
    numberstyle=\tiny\color{codegray},
    stringstyle=\color{codepurple},
    breakatwhitespace=false,
    breaklines=true,
    basewidth=0.5em,
    captionpos=b,
    columns=fixed,
    keepspaces=true,
    numbers=left,
    numbersep=5pt,
    showspaces=false,
    showstringspaces=false,
    showtabs=false,
    tabsize=4
}

\endinput
