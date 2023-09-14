
all: rlibm.a

CC = gcc

rlibmfiles = sinf cosf tanf

rlibmsrcdir = .
rlibmobjdir = obj

rlibmsrc = $(patsubst %, $(rlibmsrcdir)/%.c, $(rlibmfiles))
rlibmobj = $(patsubst %, $(rlibmobjdir)/%.o, $(rlibmfiles))

$(rlibmobjdir)/%.o: $(rlibmsrcdir)/%.c
	@mkdir -p $(rlibmobjdir)
	$(CC) -O3 -march=native  -c $^  -o $@

rlibm.a: $(rlibmobj)
	ar rc $@ $^

clean:
	rm -rf  obj *~ rlibm.a
