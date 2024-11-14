
int alectura= 34;
int umbral=(2000);

void setup() {
  // put your setup code here, to run once:
  Serial.begin (115200);

}

void loop() {
  // put your main code here, to run repeatedly:
 int lectura=analogRead(alectura);

  if (lectura>umbral)
  { 
    Serial.print(1);
  }
  else
  {
    Serial.print (0);
  }
  delay(50);
}
