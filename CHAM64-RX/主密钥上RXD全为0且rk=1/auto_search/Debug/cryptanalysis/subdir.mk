################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../cryptanalysis/DDTAB_16bit_table.c \
../cryptanalysis/DDTAB_4bit_table.c \
../cryptanalysis/DDTAB_8bit_table.c \
../cryptanalysis/HW_space_16bit_table.c \
../cryptanalysis/search.c 

OBJS += \
./cryptanalysis/DDTAB_16bit_table.o \
./cryptanalysis/DDTAB_4bit_table.o \
./cryptanalysis/DDTAB_8bit_table.o \
./cryptanalysis/HW_space_16bit_table.o \
./cryptanalysis/search.o 

C_DEPS += \
./cryptanalysis/DDTAB_16bit_table.d \
./cryptanalysis/DDTAB_4bit_table.d \
./cryptanalysis/DDTAB_8bit_table.d \
./cryptanalysis/HW_space_16bit_table.d \
./cryptanalysis/search.d 


# Each subdirectory must supply rules for building sources it contributes
cryptanalysis/%.o: ../cryptanalysis/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -fopenmp -msse4.2 -mcmodel=large -I/home/hmj110131/workshop/auto_search/variables -I"/home/hmj110131/workshop/auto_search/RX-cryptanalysis" -I/home/hmj110131/workshop/auto_search/printf -I/opt/mpich/include -I/home/hmj110131/workshop/auto_search/cryptanalysis -I/home/hmj110131/workshop/auto_search/ciphers -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


