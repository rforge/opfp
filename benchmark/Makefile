HOCKING-RIGAILL-benchmark.pdf: HOCKING-RIGAILL-benchmark.tex refs.bib figure-systemtime-arrays.pdf figure-microbenchmark-arrays.png figure-systemtime-simulation.pdf 
	pdflatex HOCKING-RIGAILL-benchmark
	bibtex HOCKING-RIGAILL-benchmark
	pdflatex HOCKING-RIGAILL-benchmark
	pdflatex HOCKING-RIGAILL-benchmark
figure-systemtime-arrays.pdf: figure-systemtime-arrays.R systemtime.arrays.RData
	R --no-save < $<
figure-microbenchmark-arrays.png: figure-microbenchmark-arrays.R microbenchmark.arrays.RData
	R --no-save < $<
figure-systemtime-simulation.pdf: figure-systemtime-simulation.R systemtime.simulation.RData
	R --no-save < $<
systemtime.arrays.RData: systemtime.arrays.R
	R --no-save < $<
microbenchmark.arrays.RData: microbenchmark.arrays.R figure-systemtime-arrays.pdf
	R --no-save < $<
systemtime.simulation.RData: systemtime.simulation.R
	R --no-save < $<