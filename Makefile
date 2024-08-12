#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = ifort
OPTSC   = -c -qopenmp -module mod -O3
OPTSL   = -qopenmp -module mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)MAIN: $(MKDIRS) $(DOBJ)main.o \
	$(DOBJ)lsupply.o \
	$(DOBJ)partest.o \
	$(DOBJ)solveactivelife.o \
	$(DOBJ)solvefirstactive.o \
	$(DOBJ)solvelastactive.o
	@rm -f $(filter-out $(DOBJ)main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MAIN

#compiling rules
$(DOBJ)bspline_kinds_module.o: src/bspline_kinds_module.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)bspline_sub_module.o: src/bspline_sub_module.f90 \
	$(DOBJ)bspline_kinds_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)glob0.o: src/glob0.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)globparams.o: src/GlobParams.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hours.o: src/hours.f90 \
	$(DOBJ)model_parameters.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)kind_module.o: src/kind_module.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)lsupply.o: src/lsupply.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)hours.o \
	$(DOBJ)root_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)main.o: src/main.f90 \
	$(DOBJ)utilities.o \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)solveinretirement.o \
	$(DOBJ)tauchen.o \
	$(DOBJ)policyfunctions_obj.o \
	$(DOBJ)simulation.o \
	$(DOBJ)statistics.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)model_parameters.o: src/Model_Parameters.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)partest.o: src/partest.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o \
	$(DOBJ)globparams.o \
	$(DOBJ)policyfunctions_obj.o \
	$(DOBJ)glob0.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)policyfunctions.o: src/PolicyFunctions.f90 \
	$(DOBJ)model_parameters.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)policyfunctions_obj.o: src/PolicyFunctions_obj.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)globparams.o \
	$(DOBJ)bspline_sub_module.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)root_module.o: src/root_module.f90 \
	$(DOBJ)bspline_kinds_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)simulation.o: src/Simulation.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)sorth.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o \
	$(DOBJ)policyfunctions_obj.o \
	$(DOBJ)globparams.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solveactivelife.o: src/SolveActiveLife.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o \
	$(DOBJ)policyfunctions_obj.o \
	$(DOBJ)globparams.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solvefirstactive.o: src/Solvefirstactive.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o \
	$(DOBJ)globparams.o \
	$(DOBJ)policyfunctions_obj.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solveinretirement.o: src/SolveInRetirement.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o \
	$(DOBJ)globparams.o \
	$(DOBJ)glob0.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solvelastactive.o: src/Solvelastactive.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o \
	$(DOBJ)globparams.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sorth.o: src/SortH.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)statistics.o: src/Statistics.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o \
	$(DOBJ)simulation.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)tauchen.o: src/tauchen.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)utilities.o: src/Utilities.f90 \
	$(DOBJ)model_parameters.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
