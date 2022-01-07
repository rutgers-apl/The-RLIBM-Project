all: floatlibm floatrnolibm 

floatlibm:floatlibm.a

floatrnolibm: floatrnolibm.a



CC = gcc

floatfiles = exp2 log2 exp log10 log cosh cosh_v2 cosh_v3 cosh_v4 sinh sinpi cospi exp10 exp10_v2
floatrnofiles = log2 log exp10 exp2 exp log10 cosh sinh cospi sinpi log-piecewise log2-piecewise exp-piecewise exp2-piecewise exp10-piecewise log10-piecewise log10-piecewise_v2


floatsrcdir = float
floatobjdir = obj

floatrnosrcdir = float-rno
floatrnoobjdir = rno-obj


floatprogsrcdir = progressive
floatprogobjdir = prog-obj

floatsrc = $(patsubst %, $(floatsrcdir)/%.c, $(floatfiles))
floatobj = $(patsubst %, $(floatobjdir)/%.o, $(floatfiles))

floatrnosrc = $(patsubst %, $(floatrnosrcdir)/%.c, $(floatrnofiles))
floatrnoobj = $(patsubst %, $(floatrnoobjdir)/%.o, $(floatrnofiles))


floatprogsrc = $(patsubst %, $(floatprogsrcdir)/%.c, $(floatprogfiles))
floatprogobj = $(patsubst %, $(floatprogobjdir)/%.o, $(floatprogfiles))


$(floatobjdir)/%.o: $(floatsrcdir)/%.c
	@mkdir -p $(floatobjdir)
	$(CC) -std=c11 -O3 -c $^ -I../include/ -o $@

$(floatrnoobjdir)/%.o: $(floatrnosrcdir)/%.c
	@mkdir -p $(floatrnoobjdir)
	$(CC) -std=c11 -O3 -c $^ -I../include/ -o $@

$(floatprogobjdir)/%.o: $(floatprogsrcdir)/%.c
	@mkdir -p $(floatprogobjdir)
	$(CC) -std=c11 -O3 -c $^ -I../include/ -o $@

floatlibm.a: $(floatobj)
	ar rc $@ $^

floatrnolibm.a: $(floatrnoobj)
	ar rc $@ $^

floatproglibm.a: $(floatprogobj)
	ar rc $@ $^

clean:
	rm -rf float/*.o float/*~ floatlibm.a floatrnolibm.a floatproglibm.a  obj *~ float-rno/*.o rno-obj float-prog/*.o prog-obj
