################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ciphers/Alzette.c \
../ciphers/CHAM.c \
../ciphers/LEA.c \
../ciphers/Siphash.c \
../ciphers/chaskey.c \
../ciphers/ciphers.c \
../ciphers/gimli.c \
../ciphers/hight.c \
../ciphers/simon.c \
../ciphers/sparx.c \
../ciphers/speck.c 

OBJS += \
./ciphers/Alzette.o \
./ciphers/CHAM.o \
./ciphers/LEA.o \
./ciphers/Siphash.o \
./ciphers/chaskey.o \
./ciphers/ciphers.o \
./ciphers/gimli.o \
./ciphers/hight.o \
./ciphers/simon.o \
./ciphers/sparx.o \
./ciphers/speck.o 

C_DEPS += \
./ciphers/Alzette.d \
./ciphers/CHAM.d \
./ciphers/LEA.d \
./ciphers/Siphash.d \
./ciphers/chaskey.d \
./ciphers/ciphers.d \
./ciphers/gimli.d \
./ciphers/hight.d \
./ciphers/simon.d \
./ciphers/sparx.d \
./ciphers/speck.d 


# Each subdirectory must supply rules for building sources it contributes
ciphers/%.o: ../ciphers/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -fopenmp -msse4.2 -mcmodel=large -I/home/hmj110131/workshop/auto_search/variables -I"/home/hmj110131/workshop/auto_search/RX-cryptanalysis" -I/home/hmj110131/workshop/auto_search/printf -I/opt/mpich/include -I/home/hmj110131/workshop/auto_search/cryptanalysis -I/home/hmj110131/workshop/auto_search/ciphers -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


