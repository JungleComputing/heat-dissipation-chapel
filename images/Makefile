# Build the images

# relies on ImageMagick's convert command
CONVERT = convert

.PHONY: 	clean all

all: pat2_8192x8192.pgm plasma_8192x8192.pgm

# Build the input images

GEN_PATTERN1 = \
	w5=`expr $$W \* 5 / 100`; \
	h5=`expr $$H \* 5 / 100`; \
	midw=`expr $$W / 2`; \
	midh=`expr $$H / 2`; \
	midw5=`expr $$midw - $$w5`; \
	midh5=`expr $$midh - $$h5`; \
	$(CONVERT) -size $${W}x$${H} xc:grey50 \
	-draw "fill grey45 rectangle 0,0 $$w5,$$H stroke grey75 line 0,0 $$W,$$H line $$midw,0 $$W,$$midh line 0,$$midh $$midw,$$H fill grey95 circle $$midw5,$$midh5 $$midw,$$midh" -compress none

GEN_PATTERN2 = \
	w5=`expr $$W \* 5 / 100`; \
	h5=`expr $$H \* 5 / 100`; \
	midw=`expr $$W / 2`; \
	midh=`expr $$H / 2`; \
	midw5=`expr $$midw - $$w5`; \
	midh5=`expr $$midh - $$h5`; \
	$(CONVERT) -size $${W}x$${H} xc:black \
	-draw "fill white circle $$midw5,$$midh5 $$midw,$$midh circle $$midw,$$h5 $$midw,0 rectangle 0,$$midh5 $$w5,$$midh" \
	-compress none

GEN_PATTERN3 = \
	w5=`expr $$W \* 5 / 100`; \
	h5=`expr $$H \* 5 / 100`; \
	midw=`expr $$W / 2`; \
	midh=`expr $$H / 2`; \
	midw5=`expr $$midw - $$w5`; \
	midh5=`expr $$midh - $$h5`; \
	$(CONVERT) -size $${W}x$${H} xc:black \
	-draw "stroke grey90 fill grey90 polygon 0,0 $$W,0 $$W,$$H fill white circle $$midw5,$$midh5 $$midw,$$midh circle $$midw,$$h5 $$midw,0 rectangle 0,$$midh5 $$w5,$$midh" \
	-compress none

areas_%.pgm:
	$(CONVERT) -size $* xc: +noise Random \
	   -virtual-pixel tile -blur 0x20 -channel G -separate -normalize +level 10%,100% \
	   -compress none $@

plasma_%.pgm:
	$(CONVERT) -size $* plasma:fractal -blur 0x10 -emboss 2 -compress none $@

gradient_%.pgm: 
	$(CONVERT) -size $* gradient:grey10-grey90 -distort SRT -45 -compress none $@

uni_%.pgm: 
	$(CONVERT) -size $* xc:grey50 -compress none $@

pat1_%.pgm: 
	W=`echo $* | cut -dx -f1`; H=`echo $* | cut -dx -f2`; $(GEN_PATTERN1) $@

pat2_%.pgm: 
	W=`echo $* | cut -dx -f1`; H=`echo $* | cut -dx -f2`; $(GEN_PATTERN2) $@

pat3_%.pgm: 
	W=`echo $* | cut -dx -f1`; H=`echo $* | cut -dx -f2`; $(GEN_PATTERN3) $@

clean: 
	rm -rf *.pgm

