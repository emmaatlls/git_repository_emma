# specify the main file
@default_files = ('master.tex');

# create a pdf version of the document with pdflatex
$pdf_mode = 1; # tex -> pdf

# specify the pdflatex options
$pdflatex = 'pdflatex -8bit -file-line-error -halt-on-error --shell-escape -synctex=1 %O %S';

# remove files with the following file extentions (as well as the default ones) when using "latexmk -c" or "latexmk -C"
$clean_ext .= "synctex.gz run.xml %R-blx.bib lol bbl";
$cleanup_includes_cusdep_generated = 1;

# Remove comment from next line to output a list of dependents after compilation ("shadows" latexmk output)
# $dependents_list = 1;

# convert svg files directly to pdf with inkscape
add_cus_dep('svg', 'pdf', 0, 'svg2pdf');
sub svg2pdf {
    $currdir = '$PWD';
    return system("inkscape --export-text-to-path --export-area-drawing --export-pdf=\"$currdir/$_[0].pdf\" \"$currdir/$_[0].svg\"");
}

# run makeindex
add_cus_dep('idx', 'ind', 0, 'makeindex');
sub makeindex{
    return system("makeindex \"$_[0].idx\"");
}
push @generated_exts, 'ind';
