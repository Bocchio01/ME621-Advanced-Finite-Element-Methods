.PHONY: r1 r2 clean

all: r2 clean

CC = xelatex

r1: ../Assignment_1/latex/Report.tex
	$(eval OUTPUT_DIR := Assignment_1)
	$(CC) -output-directory=../$(OUTPUT_DIR) -include-directory=../Assignment_1/latex/ $<
	$(CC) -output-directory=../$(OUTPUT_DIR) -include-directory=../Assignment_1/latex/ $<

r2: ../Assignment_2/latex/Report.tex
	$(eval OUTPUT_DIR := Assignment_2)
	$(CC) -output-directory=../$(OUTPUT_DIR) -include-directory=../Assignment_2/latex/ $<
	$(CC) -output-directory=../$(OUTPUT_DIR) -include-directory=../Assignment_2/latex/ $<


clean:
	del /Q ..\$(OUTPUT_DIR)\*.aux
	del /Q ..\$(OUTPUT_DIR)\*.lof
	del /Q ..\$(OUTPUT_DIR)\*.log
	del /Q ..\$(OUTPUT_DIR)\*.lot
	del /Q ..\$(OUTPUT_DIR)\*.out
	del /Q ..\$(OUTPUT_DIR)\*.toc
	del /Q ..\$(OUTPUT_DIR)\*.fdb_latexmk
	del /Q ..\$(OUTPUT_DIR)\*.fls
	del /Q ..\$(OUTPUT_DIR)\*.gz