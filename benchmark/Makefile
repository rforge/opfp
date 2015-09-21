Benchmark_ForPaper.pdf: Benchmark_ForPaper.tex refs.bib figure-systemtime-simulationLarge.pdf figure-systemtime-simulation.pdf figure-systemtime-arrays.pdf table-neuroblastoma-test-error.tex
	pdflatex Benchmark_ForPaper
	bibtex Benchmark_ForPaper
	pdflatex Benchmark_ForPaper
	pdflatex Benchmark_ForPaper
table-neuroblastoma-test-error.tex: table-neuroblastoma-test-error.R algo.colors.R
	R --no-save < $<
figure-systemtime-simulationLarge.pdf: figure-systemtime-simulationLarge.R systemtime.simulationLarge.RData algo.colors.R
	R --no-save < $<
systemtime.simulationLarge.RData: systemtime.simulationLarge.R
	R --no-save < $<
HOCKING-RIGAILL-benchmark.pdf: HOCKING-RIGAILL-benchmark.tex refs.bib figure-systemtime-arrays.pdf figure-microbenchmark-arrays.png figure-systemtime-simulation.pdf 
	pdflatex HOCKING-RIGAILL-benchmark
	bibtex HOCKING-RIGAILL-benchmark
	pdflatex HOCKING-RIGAILL-benchmark
	pdflatex HOCKING-RIGAILL-benchmark
figure-systemtime-arrays.pdf: figure-systemtime-arrays.R systemtime.arrays.RData algo.colors.R
	R --no-save < $<
figure-microbenchmark-arrays.png: figure-microbenchmark-arrays.R microbenchmark.arrays.RData
	R --no-save < $<
figure-systemtime-simulation.pdf: figure-systemtime-simulation.R systemtime.simulation.RData algo.colors.R
	R --no-save < $<
systemtime.arrays.RData: systemtime.arrays.R
	R --no-save < $<
microbenchmark.arrays.RData: microbenchmark.arrays.R 
	R --no-save < $<
systemtime.simulation.RData: systemtime.simulation.R
	R --no-save < $<