################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../RX-cryptanalysis/RX_Alzette.c \
../RX-cryptanalysis/RX_CHAM.c \
../RX-cryptanalysis/RX_Chaskey.c \
../RX-cryptanalysis/RX_Siphash.c \
../RX-cryptanalysis/RX_crypt.c \
../RX-cryptanalysis/RX_speck_data.c \
../RX-cryptanalysis/RX_speck_key.c \
../RX-cryptanalysis/XOR_offset_table.c 

OBJS += \
./RX-cryptanalysis/RX_Alzette.o \
./RX-cryptanalysis/RX_CHAM.o \
./RX-cryptanalysis/RX_Chaskey.o \
./RX-cryptanalysis/RX_Siphash.o \
./RX-cryptanalysis/RX_crypt.o \
./RX-cryptanalysis/RX_speck_data.o \
./RX-cryptanalysis/RX_speck_key.o \
./RX-cryptanalysis/XOR_offset_table.o 

C_DEPS += \
./RX-cryptanalysis/RX_Alzette.d \
./RX-cryptanalysis/RX_CHAM.d \
./RX-cryptanalysis/RX_Chaskey.d \
./RX-cryptanalysis/RX_Siphash.d \
./RX-cryptanalysis/RX_crypt.d \
./RX-cryptanalysis/RX_speck_data.d \
./RX-cryptanalysis/RX_speck_key.d \
./RX-cryptanalysis/XOR_offset_table.d 


# Each subdirectory must supply rules for building sources it contributes
RX-cryptanalysis/%.o: ../RX-cryptanalysis/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -fopenmp -msse4.2 -mcmodel=large -I/home/hmj110131/workshop/auto_search/variables -I"/home/hmj110131/workshop/auto_search/RX-cryptanalysis" -I/home/hmj110131/workshop/auto_search/printf -I/opt/mpich/include -I/home/hmj110131/workshop/auto_search/cryptanalysis -I/home/hmj110131/workshop/auto_search/ciphers -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


