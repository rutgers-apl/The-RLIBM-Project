all: rlibm.a


CC = gcc

rlibmfiles = log2 log exp10 exp2 exp log10 cosh sinh cospi sinpi 



rlibmsrcdir = .
rlibmobjdir = obj


rlibmsrc = $(patsubst %, $(rlibmsrcdir)/%.c, $(rlibmfiles))
rlibmobj = $(patsubst %, $(rlibmobjdir)/%.o, $(rlibmfiles))



$(rlibmobjdir)/%.o: $(rlibmsrcdir)/%.c
	@mkdir -p $(rlibmobjdir)
	$(CC) -std=c11 -O3 -c $^  -o $@


rlibm.a: $(rlibmobj)
	ar rc $@ $^


clean:
	rm -rf  obj *~ rlibm.a
