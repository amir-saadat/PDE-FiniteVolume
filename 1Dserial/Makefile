#------------------------ Makefile for serial code --------------------#
##usage:  make compile ; make run  
#----------------------------------------------------------------------#
#If your compiler is NOT on your path (for your shell) then 
#   you need to insert the full path, e.g. /opt/intel/..../bin/ifort
#>>>>>>>>>>>>>> set appropriate COMPiler:  uncommend one
 COMP   = gcc
 FLAGs  = -g
 #FLAGs = -O3 
 LFLAGs = -lm
 EXT    = c

# COMP   = gfortran
# FLAGs  = -g
##FLAGs  = -O3 
# LFLAGs = 
# EXT    = f
#========================== set source code  ==========================#
#>>>>>>>>>>>>>>>>>>>> set names for PROG and std I/O : 
 PROG = code1D
#PROG = code2D 
 INPUT    = data.dat
 OUTPUT   = o.out
#>>>>>>>>>>>>>>>>>>>> set code components: 
  CODE_o = z.main.o  z.io.o  z.setup.o  z.update.o 

##-------------------- not needed for most 'make', does no harm:
# CODE   = z.main.$(EXT)  z.io.$(EXT)  z.setup.$(EXT)  z.update.$(EXT)
 
#======================= create executable: make compile ==============#
#------------> lines below a directive MUST start with TAB <-----------#
$(CODE_o):%.o: %.$(EXT)
	$(COMP) $(FLAGs) -c $< -o $@
 
compile: $(CODE_o) 
	@echo "         compiling with  $(COMP) :" 
	$(COMP) $(FLAGs)  $(CODE_o)  -o $(PROG).x  $(LFLAGs)
	@echo " >>> compiled on `hostname -s` with  $(COMP) <<<" 
	rm -f z.*.o 
#============================ execute: make run ==========# 
run: 
	./$(PROG).x < $(INPUT) > $(OUTPUT) 
#====================================================================# 
clean: 
	rm -f  o.*  DONE 
#---------------------------------------------------------------------- 

