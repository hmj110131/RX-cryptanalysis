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
	gcc -I"/home/hmj110131/workshop/auto_search/ciphers" -I/home/hmj110131/workshop/auto_search/printf -I/home/hmj110131/workshop/auto_search/variables -I/home/hmj110131/workshop/auto_search/cryptanalysis -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


