.PHONY: reports r1 clean

all: reports clean

CC = xelatex

r1: ../Assignment_1/latex/Report.tex
	$(CC) -output-directory=../Assignment_1/latex -include-directory=../Assignment_1/latex/ $<

reports: r1 # r2 r3 r4

clean:
	del /Q *.aux
	del /Q *.lof
	del /Q *.log
	del /Q *.lot
	del /Q *.out
	del /Q *.toc
	del /Q *.fdb_latexmk
	del /Q *.fls
	del /Q *.gz