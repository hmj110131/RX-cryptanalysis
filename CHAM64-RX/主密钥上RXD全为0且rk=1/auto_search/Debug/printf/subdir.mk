################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../printf/printinfo.c 

OBJS += \
./printf/printinfo.o 

C_DEPS += \
./printf/printinfo.d 


# Each subdirectory must supply rules for building sources it contributes
printf/%.o: ../printf/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -fopenmp -msse4.2 -mcmodel=large -I/home/hmj110131/workshop/auto_search/variables -I"/home/hmj110131/workshop/auto_search/RX-cryptanalysis" -I/home/hmj110131/workshop/auto_search/printf -I/opt/mpich/include -I/home/hmj110131/workshop/auto_search/cryptanalysis -I/home/hmj110131/workshop/auto_search/ciphers -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


