################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include variables/subdir.mk
-include printf/subdir.mk
-include cryptanalysis/subdir.mk
-include ciphers/subdir.mk
-include RX-cryptanalysis/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: auto_search

# Tool invocations
auto_search: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc -fopenmp -msse4.2 -mcmodel=large -L/opt/mpich/lib -o "auto_search" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) auto_search
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets
