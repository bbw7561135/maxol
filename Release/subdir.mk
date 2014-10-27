################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bound_cond.cpp \
../evolution_electric.cpp \
../evolution_magnetic.cpp \
../init_cond.cpp \
../maxol.cpp \
../output.cpp 

OBJS += \
./bound_cond.o \
./evolution_electric.o \
./evolution_magnetic.o \
./init_cond.o \
./maxol.o \
./output.o 

CPP_DEPS += \
./bound_cond.d \
./evolution_electric.d \
./evolution_magnetic.d \
./init_cond.d \
./maxol.d \
./output.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


