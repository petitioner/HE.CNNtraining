################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CNNtraining.cpp \
../src/Methods.cpp \
../src/Tools.cpp 

OBJS += \
./src/CNNtraining.o \
./src/Methods.o \
./src/Tools.o 

CPP_DEPS += \
./src/CNNtraining.d \
./src/Methods.d \
./src/Tools.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/home/john/eclipse-workspace/MyHECNNtraining/HEAAN/HEAAN/src" -I"/home/john/eclipse-workspace/MyHECNNtraining/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


