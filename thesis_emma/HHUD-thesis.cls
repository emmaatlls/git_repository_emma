\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{hhud-thesis}[2024/12/09]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrbook}}

\ProcessOptions
\LoadClass[
a4paper,%
12pt,%
headsepline,%
chapterprefix,%
parskip,%
listof=totoc,%
index=totoc,%
bibliography=totoc,%
cleardoublepage=empty,%
DIV=14,%
]
{scrbook}

\PassOptionsToPackage{hyphens}{url} % breaks long URLs at hyphens % Do this early, so that if a package loads 'url', it will use this option.

\RequirePackage[T1]{fontenc} % Gibt Unicode Zeichen korrekt aus. (ö anstelle von "o)
\RequirePackage[utf8]{inputenc} % erlaubt direkte Eingabe von Unicode-Zeichen im Latex Code, insbesondere Umlaute % technisch nicht mehr notwendig, aber für die Sicherheit
\RequirePackage{babel} % Ersetzt "Section" mit "Abschnitt" etc. je nach Sprache
\RequirePackage[babel,tracking]{microtype}

%% Fonts 
\RequirePackage{newtxtext} % Times New Roman
\RequirePackage{newtxmath} % Times New Roman Math % Also loads amsmath and amssymb for extended math support
\RequirePackage[scaled=0.92]{helvet} % Set sans-serif font to Helvetica
\RequirePackage{courier} % Set monospaced font to Courier

\RequirePackage{ifpdf}
\RequirePackage{color}
\RequirePackage{setspace}
\RequirePackage[plainpages=false,pdfusetitle,pdfpagelabels]{hyperref}
\ifpdf
\definecolor{brown}{cmyk}{0, 0.81, 1, 0.60}
\hypersetup{
	colorlinks=true,
	linkcolor=black,
	citecolor=black,
	urlcolor=black,
	bookmarksnumbered=true,
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage[style=apa, backend=biber]{biblatex}

%%%% Workarounds

\AddToHook{begindocument/before}{\@ifpackageloaded{minted}{\removefromtoclist[float]{lol}}{}}

\def\displaymath{\typeout{*** Do not use DISPLAYMATH. Use \[...\] instead ***}}
\def\eqnarray{\typeout{*** Do not use EQNARRAY(*). Use align(*)-environment instead ***}}


%%%%
\newcommand{\thesistypegerman}{\graduationlevel{}arbeit}
%% Title Page
\renewcommand{\maketitle}{%
	\pdfbookmark[0]{\IfLanguageName{ngerman}{Titelseite}{Front Page}}{frontpage}
	
\begin{titlepage}

	\noindent
	\begin{minipage}[t]{0.48\textwidth}
		\raggedright
		\IfLanguageName{ngerman}
		{\includegraphics[width=\linewidth]{static/HHU_Logo_Wort-Bild-Marke_horizontal_RGB}}
		{\includegraphics[width=\linewidth]{static/Logo_HHU_+Name_EN_horizontal_rgb_+Safezone}}
	\end{minipage}%
	\hfill
	\begin{minipage}[t]{0.48\textwidth}
		\raggedleft
		\IfLanguageName{ngerman}
		{\includegraphics[width=\linewidth]{static/HHU_Logo_Wort-Bild-Marke_horizontal_RGB}}
		{\includegraphics[width=\linewidth]{static/UniversitaetzuKoelnLogoENG.jpg}}
	\end{minipage}
	
	\vspace{45mm}
		
		% Titel
		\begin{center}
			\textbf{\huge \@title}
			\\
			
			
			\vspace{10mm}
			
			\textbf{\Large \@author}
			
			\begin{tabular}[t]{ll}
				\IfLanguageName{ngerman}{Matrikelnummer}{Student ID}: & \studentid \\
			\end{tabular}
		\end{center}
		
		\vspace{50mm}
		
		\begin{center}
			
			{\Large \thesistype}
			
			\IfLanguageName{ngerman}{Vorgelegt am Lehrstuhl für Rechnernetze.}{Submitted to the Institute of Biology at Heinrich-Heine-University, Düsseldorf.}
			
			\vfill
			
			\begin{tabular}[t]{ll}
				\IfLanguageName{ngerman}{Abgabe der Arbeit}{Submission}: & \DTMmonthname{\thesissubmissionmonth} \thesissubmissionyear \\
				\ifthenelse{\NOT \equal{\graduationlevel}{Seminar}}
				{% Abschlussarbeit (Print reviewers for Bachelor & Master)
					\IfLanguageName{ngerman}{Erstgutachten}{First Reviewer}:         & \firstreviewer \\
					\IfLanguageName{ngerman}{Zweitgutachten}{Second Reviewer}:       & \secondreviewer \\
				}{}
			\end{tabular}
		\end{center}
	\end{titlepage}
}



\DefineBibliographyStrings{ngerman}{%
	mathesis = {Masterarbeit},
	bathesis = {Bachelorarbeit},
	phdthesis = {Dissertation}
}

\DefineBibliographyStrings{english}{%
	mathesis = {Master's Thesis},
	bathesis = {Bachelor's Thesis},
	phdthesis = {PhD Thesis}
}

\fi