#include <avr/io.h>
#include <avr/interrupt.h>
//Using the  AVR libraries. Compiled in C these unlock the potential of the AVR microcontrollers upon which Arduino is based.

#define Divisions (200)// Sub divisions of sisusoidal wave.
static unsigned int lookUp_A[Divisions];
static unsigned int lookUp_B[Divisions];
static unsigned int lookUp_C[Divisions];
//Preparing each look up table as an array with "Divisions" number of entries.

//Defining values that will be used later

static int MHz = 16;                                   // Micro clock frequency
static int freq = 10;                                        // Sinusoidal frequency
uint32_t DutyCycle;                                         // The length of each PWM Duty Cycle
uint32_t period;                                            // Period of each sinosoid waveform
uint32_t halfperiod;                                        // Half of the period of each sinosoid waveform
uint32_t Third;                                             //Also known as 1/3 AKA 0.3333333333333333333333
uint32_t TwoThirds;                                         //Also known as 2/3 AKA 0.666666666666666666667
uint32_t X;                                                 //Used as a counting variable in the loop for Phase A
uint32_t Y;                                                 //Used as a counting variable in the loop for Phase B
uint32_t Z;                                                 //Used as a counting variable in the loop for Phase C
float ClocksperSecond = 62500;                              // Including pre-scaler, how many system clocks will occur each second
float ClocksperPeriod = ClocksperSecond / freq;             // How many clocks will happen per period of sinosoid waveforn
float ClocksperCycle = 65536;                               //System clocks per system cycle. Will only change with the pre-scaler.
float PeriodClockFraction = ClocksperPeriod / ClocksperCycle; //What is the fraction of the duty cycle against length of the system cycle? Measured in system clocks.
const double currentFactor = 2 / 3.3;                       // 3.3V per 2A which is 0.909

// Uno Pins
//If you want to run single-phase AC on an Arduino Uno (and Arduino Motor Shield Rev 3)

//[Physical Pin Number] & Pin as denoted on Arduino IDE
//const int MotorPinA = PB4;                                //Pin [18] 12
//const int MotorSpeedPinA = PD3;                           //Pin [5] 3
//const int MotorBrakePinA = PB1;                           //Pin [15] 9
//const int CurrentSensePinA = A0;                           //Current Sensing

//Mega pins
//If you want to run three-phase AC on an Arduino Mega (and three Arduino Unos and three Arduino Motor Shields Rev 3's )

//[Physical Pin Number] & Pin as denoted on Arduino IDE
const int MotorPinA = PH1;                                  //Pin [13] & 16 MotorPinA for MotorShield switches the direction of the current and generates the negative component of the waveform
const int MotorSpeedPinA = PH6;                             //Pin [18] & 9 MotorSpeedPinA for MotorShield provides a PWM signal, and will allow modulation of the duty cycle
const int MotorBrakePinA = PH5;                             //Pin [17] & 8 MotorBrakePinA for MotorShield provides will cut the signal if necessary (i.e. Circuit breaker with CurrentSensing)
const int CurrentSensePinA = A0;                            //Current Sensing

const int MotorPinB = PA5;                                  //Pin [73] & 27
const int MotorSpeedPinB = PB7;                             //Pin [26] & 13
const int MotorBrakePinB = PA0;                             //Pin [78] & 22
const int CurrentSensePinB = A1;                            //Current Sensing

const int MotorPinC      = PL7;                             //Pin [42] & 42
const int MotorSpeedPinC = PB6;                             //Pin [25] & 12
const int MotorBrakePinC = PL6;                             //Pin [41] & 43
const int CurrentSensePinC = A2;                            //Current Sensing


uint32_t compA       =  PeriodClockFraction * ClocksperCycle / 2;                         //Interrupt compare value A, (halfway through cycle)
uint32_t compB       =  PeriodClockFraction * ClocksperCycle;                             //Interrupt compare value B, at the end of the cycle
uint32_t TimerBStart =  PeriodClockFraction * 0 * ClocksperCycle;                  //Offset for Phase A
uint32_t TimerCStart =  PeriodClockFraction * ClocksperCycle * 1 / 3;              //Offset for Phase B
uint32_t TimerAStart =  PeriodClockFraction * ClocksperCycle * 2 / 3;              //Offset for Phase C


// Counter and Compare Values
//For Uno
//uint32_t t1_load = TimerAStart;
//uint32_t t1_compA = compA;
//uint32_t t1_compB = compB;

//For Mega
uint32_t t3_load = TimerAStart;
uint32_t t3_compA = compA;
uint32_t t3_compB = compB;

uint32_t t4_load = TimerBStart;
uint32_t t4_compA = compA;
uint32_t t4_compB = compB;

uint32_t t5_load = TimerCStart;
uint32_t t5_compA = compA;
uint32_t t5_compB = compB;           //seems you can only really trust a uint32_t, dunno why int const in and float didnt work :/

void setup() {
  cli();                             //stop interrupts while we set up the timer
  Serial.begin(9600);                //seial monitor initialized, use Serial.print to see what is happening to variables within the sketch.

  double temp;                       //Double varible for <math.h> functions.
  period = ClocksperPeriod;
  halfperiod = period / 2;
  DutyCycle = period / Divisions;    //Period in Microseconds, remember it has a resolution of 4.


  //Lookup Table for Phase 1/A
  for (int i = 0; i < Divisions; i++) {                   // Filling each value in the array as its corresponding index. i.e. [1] = 1 [2] = 2 [3] = 3 etc.
    temp = abs(sin((i * 2 * M_PI / (Divisions)))) * 255;  //Values in array now correspond to a sine wave where the peak values is now 255. 255 later corresponding to 100% duty cycle
    lookUp_A[i] = (int)(temp + 0.5);                      // Round up to integer.
    lookUp_B[i] = (int)(temp + 0.5);  
    lookUp_C[i] = (int)(temp + 0.5);  
    Serial.print(i);
    Serial.print("  :  ");
        Serial.println(lookUp_A[i]);
  }
  
  //DDRB is a bitmap oriented register that controls the Data Direction of each bit on Port B. i.e. 0 for Read and 1 for Write

  //Uno
  //  DDRB |= (1 << MotorPinA);
  //  DDRD |= (1 << MotorSpeedPinA);
  //  DDRB |= (1 << MotorBrakePinA);

  //Mega
  DDRH |= (1 << MotorPinA);
  DDRH |= (1 << MotorSpeedPinA);
  DDRH |= (1 << MotorBrakePinA);

  DDRA |= (1 << MotorPinB);
  DDRB |= (1 << MotorSpeedPinB);
  DDRA |= (1 << MotorBrakePinB);

  DDRL |= (1 << MotorPinC);
  DDRB |= (1 << MotorSpeedPinC);
  DDRL |= (1 << MotorBrakePinC);

  //PORTB is the register the code uses to set the port pins of Port B if writing, i.e. 0 for LOW and 1 for HIGH

  //Uno
  //  PORTB |= (1 << MotorPinA);
  //  PORTD |= (1 << MotorSpeedPinA);
  //  PORTB |= (0 << MotorBrakePinA);

  //Mega
  PORTH |= (1 << MotorPinA);
  PORTH |= (1 << MotorSpeedPinA); //This acts like analogWrite(255) i.e. 100% duty cycle
  PORTH |= (0 << MotorBrakePinA);

  PORTA |= (1 << MotorPinB);
  PORTB |= (1 << MotorSpeedPinB);
  PORTA |= (0 << MotorBrakePinB);

  PORTL |= (1 << MotorPinC);
  PORTB |= (1 << MotorSpeedPinC);
  PORTL |= (0 << MotorBrakePinC);

  //Uno
  //TIFR1 |= (1 << OCF1B);
  //TIFR1 |= (1 << OCF1A);

  //   TCCR1A = 0;

  //Mega
  TIFR3 |= (1 << OCF3B);
  TIFR3 |= (1 << OCF3A);
  //Enabling Timer Flags for Outut compare A and B. i.e. the type of interrupt we will be using.

  TCCR3A = 0;

  // Clearing TCCR3A
  //Bits in TCCR3A can be set to do a bunch of things, but we don't need any of them for this.
  //
  TIFR4 |= (1 << OCF4B);
  TIFR4 |= (1 << OCF4A);

  // Clearing TCCR4A
  TCCR4A = 0;

  TIFR5 |= (1 << OCF5B);
  TIFR5 |= (1 << OCF5A);

  // Clearing TCCR5A
  TCCR5A = 0;
  // Reset Timer1 Control Reg B
  //  TCCR1B = 0;
  //   Reset Timer3 Control Reg B
  TCCR3B = 0;
  //  Reset Timer4 Control Reg B
  TCCR4B = 0;
  // Reset Timer5 Control Reg B
  TCCR5B = 0;
  //Just making sure this zero before we set any bits

  //set pre-scalar to 256
  //Arduino clock operates at 16,000,000 Hz. 16 bit timers roll over after 65536 clocks (2^16). In order to use frequencies like 1Hz or 50Hz a prescaler is applied.
  //This means that 256 system clocks will pass everytime a single clock is counted by any of the timers

  //Timer1B

  //  TCCR1B |= (1 << CS12);
  //  TCCR1B |= (0 << CS11);
  //  TCCR1B |= (0 << CS10);

  //Timer3B

  TCCR3B |= (1 << CS32);
  TCCR3B |= (0 << CS31);
  TCCR3B |= (0 << CS30);
  //
  //   //Timer4B
  TCCR4B |= (1 << CS42);
  TCCR4B |= (0 << CS41);
  TCCR4B |= (0 << CS40);
  //
  //   //Timer5B
  TCCR5B |= (1 << CS52);
  TCCR5B |= (0 << CS51);
  TCCR5B |= (0 << CS50);
  //Setting these bits correspond to pre-scaler of 256


  //Enable the Timer interrupts to take place
  //TCNT1 is Timer1
  //TCNT1 = t1_load;
  //OCR1A = t1_compA; //Comp A
  //OCR1B = t1_compB; //Comp B

  //  //TCNT3 is Timer3
  TCNT3 = t3_load;
  OCR3A = t3_compA; //Comp A
  OCR3B = t3_compB; //Comp B
  //
  ////    //TCNT4 is Timer4
  TCNT4 = t4_load;
  OCR4A = t4_compA; //Comp A
  OCR4B = t4_compB; //Comp B
  ////
  ////    //TCNT5 is Timer5
  TCNT5 = t5_load;
  OCR5A = t5_compA; //Comp A
  OCR5B = t5_compB; //Comp B

  sei(); //enable global interrupts
  //  //Timer1
  //    TIMSK1 = 0;
  //  TIMSK1 |= (1 << OCIE1A);
  //  TIMSK1 |= (1 << OCIE1B);
  //  //Timer3
  TIMSK3 = 0;
  TIMSK3 |= (1 << OCIE3A);
  TIMSK3 |= (1 << OCIE3B);
  ////Specifically allows Timer 3 to intitate interrupts
  //// //Timer4
  TIMSK4 = 0;
  TIMSK4 |= (1 << OCIE4A);
  TIMSK4 |= (1 << OCIE4B);
  ////Specifically allows Timer 4 to intitate interrupts
  ////Timer5
  TIMSK5 = 0;
  TIMSK5 |= (1 << OCIE5A);
  TIMSK5 |= (1 << OCIE5B);  //Disable or enable sets of timers using this timsk stuff
  //Specifically allows Timer 4 to intitate interrupts
}

void loop() {
//  double currentSenseA = map(analogRead(CurrentSensePinA), 0, 1024, 0, 5000);
//  double currentSenseB = map(analogRead(CurrentSensePinB), 0, 1024, 0, 5000);
//  double currentSenseC = map(analogRead(CurrentSensePinC), 0, 1024, 0, 5000);
//  double CurrentA = currentSenseA * currentFactor;                                         // get channel A current from voltage and current factor (from data sheet)
//  double CurrentB = currentSenseB * currentFactor;                                         // get channel A current from voltage and current factor (from data sheet)
//  double CurrentC = currentSenseC * currentFactor;
//  // get channel A current from voltage and current factor (from data sheet)
//  if (CurrentA > 2000
//      || CurrentB > 2000
//      || CurrentC > 2000
//     )
//  {
//    //    PORTB |= (1 << MotorBrakePinA);
//    PORTH |= (1 << MotorBrakePinA);
//    PORTA |= (1 << MotorBrakePinB);
//    PORTL |= (1 << MotorBrakePinC);
//  }

  //       X = int((TCNT1/DutyCycle)+0.5);                               //X represents the number of DutyCycles that have passed on the current Timer 3 clock cycle
  //     analogWrite(3, lookUp_A[X]);                                  //Assigns the corresponding sine value to analogwrite, changing the duty cycle as appropriate
  //Serial.print(lookUp[X]);                                               //Records the clock counter of the timer as this process occurs.
  //Serial.print("---");

  X = int(((TCNT3) / DutyCycle));                           //X represents the number of DutyCycles that have passed on the current Timer 3 clock cycle
  analogWrite(9, lookUp_A[X]);                                  //Assigns the corresponding sine value to analogwrite, changing the duty cycle as appropriate
//  Serial.print("Timer 3: ");
//  Serial.print(TCNT3);                                               //Records the clock counter of the timer as this process occurs.
//  Serial.print("---");
//  Serial.println(lookUp[X]);

  Y = int(((TCNT4)/ DutyCycle));
  analogWrite(13, lookUp_B[Y]);
//  Serial.print("Timer 4: ");
//  Serial.print(TCNT4);
//  Serial.print("---");
//   Serial.println(lookUp[Y]);

  Z = int(((TCNT5) / DutyCycle));
  analogWrite(12, lookUp_C[Z]);
//  Serial.print("Timer 5: ");
//  Serial.print(TCNT5);
//    Serial.print("---");
//   Serial.println(lookUp[Z]);
} //end loop


//Interrupts - These are called when the interrupt conditions establish previously are met
//Uno
//ISR(TIMER1_COMPA_vect){
//    PORTB ^= (1 << MotorPinA); //Current direction changed
////  PORTD ^= (1 << MotorSpeedPinA);
////  PORTB ^= (1 << MotorBrakePinA);
////  Serial.print(TCNT1);
////  Serial.println(lookUp[x]);
//}

//ISR(TIMER1_COMPB_vect){
//    PORTB ^= (1 << MotorPinA); //Current direction changed
////  PORTD ^= (1 << MotorSpeedPinA);
////  PORTB ^= (1 << MotorBrakePinA);
////  Serial.print(TCNT1);
////  Serial.println(lookUp[X]);
//    TCNT1 = 0;
//}

//Mega
ISR(TIMER3_COMPA_vect) {
  PORTH ^= (1 << MotorPinA); //Current direction changed
  //  PORTH ^= (1 << MotorSpeedPinA);
  //  PORTH ^= (1 << MotorBrakePinA);
    Serial.print("Interrupt 3A: ");
    Serial.print(TCNT3);
    Serial.print("  :  ");
  Serial.print(X);
  Serial.print("  :  ");
    Serial.println(lookUp_A[X]);
}
//
ISR(TIMER3_COMPB_vect) { //Current direction changed
  PORTH ^= (1 << MotorPinA);
//    PORTH ^= (1 << MotorSpeedPinA);
//    PORTH ^= (1 << MotorBrakePinA);
    Serial.print("Interrupt 3B: ");
     Serial.print(TCNT3);
       Serial.print("  :  ");
      Serial.print(X);
      Serial.print("  :  ");
      Serial.println(lookUp_A[X]);
  TCNT3 = 0;
}
////
ISR(TIMER4_COMPA_vect) { //Current direction changed
  PORTA ^= (1 << MotorPinB);
  //PORTLB ^= (1 << MotorSpeedPinB);
  // PORTA ^= (1 << MotorBrakePinB);
    Serial.print("Interrupt 4A: ");
     Serial.print(TCNT4);
       Serial.print("  :  ");
    Serial.print(Y);
    Serial.print("  :  ");
      Serial.println(lookUp_B[Y]);
}
ISR(TIMER4_COMPB_vect) { //Current direction changed
  PORTA ^= (1 << MotorPinB);
  //  PORTB ^= (1 << MotorSpeedPinB);
  //  PORTA ^= (1 << MotorBrakePinB);
    Serial.print("Interupt 4B: ");
     Serial.print(TCNT4);
       Serial.print("  :  ");
   Serial.print(Y);
   Serial.print("  :  ");
      Serial.println(lookUp_B[Y]);
  TCNT4 = 0;
}
////
ISR(TIMER5_COMPA_vect) { //Current direction changed
  PORTL ^= (1 << MotorPinC);
  //  PORTB ^= (1 << MotorSpeedPinC);
  //  PORTL ^= (1 << MotorBrakePinC);
    Serial.print("Interupt 5A: ");
     Serial.print(TCNT5);
       Serial.print("  :  ");
    Serial.print(Z);
           Serial.print("  :  ");
      Serial.println(lookUp_C[Z]);
}

ISR(TIMER5_COMPB_vect) { //Current direction changed
  PORTL ^= (1 << MotorPinC);
  //  PORTB ^= (1 << MotorSpeedPinC);
  //  PORTL ^= (1 << MotorBrakePinC);
  Serial.print("Interupt 5B: ");
      Serial.print(TCNT5);
       Serial.print("  :  ");
      Serial.print(Z);
       Serial.print("  :  ");
      Serial.println(lookUp_C[Z]);
  TCNT5 = 0;
}
