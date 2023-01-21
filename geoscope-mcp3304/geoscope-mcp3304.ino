// Geoscope, V. Ziemann, 230121

#include "fft.c"
#define npts 512
int bits;
complx wave[npts]; 

const char* ssid     = "YOUR_SSID";
const char* password = "YOUR_PASSWORD";
const int port=81;

uint8_t channel=0;

#include <ESP8266WiFi.h>
#include <ESP8266WebServer.h>
#include <WebSocketsServer.h>
#include <ArduinoJson.h>
#include "FS.h"
ESP8266WebServer server2(80); 
WebSocketsServer webSocket = WebSocketsServer(port);

#include <Ticker.h>
Ticker SampleFast, SampleSlow;
volatile uint16_t isamp=0;
int sample_period=10,sample_buffer_ready=0,spectrum=1,vscale=1;
unsigned long time_to_wait=50;
int sample_buffer[npts];
volatile uint8_t websock_num=0;
char out[4000]; DynamicJsonDocument doc(4000);

void sampleslow_action(uint8_t num) {
  SampleSlow.detach();
  SampleFast.attach_ms(sample_period,samplefast_action,num);
}

void samplefast_action(uint8_t num) {
//sample_buffer[isamp]=analogRead(0);
  sample_buffer[isamp]=read_adc(channel);
  isamp++;
  if (npts == isamp) {   
     SampleFast.detach();
     isamp=0;
     sample_buffer_ready=1;
  }
}

//...................................................SPI ADC
#include <SPI.h>
#define CS 15
int read_adc(int channel) {  //...0 to 3.......read_adc
  int adcvalue=0, b1=0, hi=0, lo=0, sign=0, reading;
  digitalWrite (CS, LOW);
  byte commandbits = B00001000; // Startbit+(diff=0)
  commandbits |= channel & 0x03;
  SPI.transfer(commandbits);
  b1 = SPI.transfer(0x00);     // always D0=0
  sign = b1 & B00010000;
  hi = b1 & B00001111;
  lo = SPI.transfer(0x00);     // input is don't care
  digitalWrite(CS, HIGH);
  reading = (hi << 8) + lo;
  if (sign) {reading = (4096 - reading) * -1;  }
  return reading;
}

//...................................................................
void webSocketEvent(uint8_t num, WStype_t type, uint8_t * payload, size_t length) {
  Serial.printf("webSocketEvent(%d, %d, ...)\r\n", num, type);
  websock_num=num;
  switch(type) {
    case WStype_DISCONNECTED:
      Serial.printf("[%u] Disconnected!\r\n", num);
      break;
    case WStype_CONNECTED: 
      {
        IPAddress ip = webSocket.remoteIP(num);
        Serial.printf("[%u] Connected from %d.%d.%d.%d url: %s\r\n", 
          num, ip[0], ip[1], ip[2], ip[3], payload);
      }
      break;
    case WStype_TEXT:
    {
      Serial.printf("[%u] get Text: %s\r\n", num, payload);
      //........................................parse JSON
      DynamicJsonDocument root(300);
      deserializeJson(root,payload);
      const char *cmd = root["cmd"];
      const long val = root["val"];
//    Serial.println(cmd); Serial.println(val);
      if (strstr(cmd,"FASTON")) { 
        sample_period=val;
        SampleFast.attach_ms(sample_period,samplefast_action,num);
      } else if (strstr(cmd,"FASTOFF")) { 
        SampleFast.detach();
        SampleSlow.detach();
      } else if (strstr(cmd,"SPECTRUM")) { 
        Serial.println("Select spectrum or time domain display");
        spectrum=val;
      } else if (strstr(cmd,"SCALE")) {
        vscale=val;
      } else if (strstr(cmd,"TTW")) {
        time_to_wait=val;
      } else if (strstr(cmd,"CHANNEL")) {
        channel=val;
      } else {
        Serial.println("Unknown command");
      }
      break;
    }
    case WStype_BIN:
      Serial.printf("[%u] get binary length: %u\r\n", num, length);
      hexdump(payload, length);
      webSocket.sendBIN(num, payload, length);  // echo back
      break;
    default:
      Serial.printf("Invalid WStype [%d]\r\n", type);
      break;
  }
}
//..................................................................setup
void setup() {
  int i;
  pinMode(CS,OUTPUT);
  digitalWrite(CS,HIGH); 
  SPI.begin();
  SPI.setFrequency(2100000);
  SPI.setBitOrder(MSBFIRST);
  SPI.setDataMode(SPI_MODE0); 
  pinMode(LED_BUILTIN,OUTPUT);
  digitalWrite(LED_BUILTIN,LOW);
  Serial.begin(115200);
  delay(1000);
  WiFi.begin(ssid, password);
  while(WiFi.status() != WL_CONNECTED) {Serial.print("."); delay(500);}
  Serial.print("\nConnected to ");  Serial.print(ssid);
  Serial.print(" with IP address: "); Serial.println(WiFi.localIP());
  webSocket.begin();
  webSocket.onEvent(webSocketEvent);
  if (!SPIFFS.begin()) {
    Serial.println("ERROR: no SPIFFS filesystem found");
  } else {
    server2.begin();          
    server2.serveStatic("/", SPIFFS, "/webpage.html");
    Serial.println("SPIFFS file system found and server started on port 80");   
  }
  digitalWrite(LED_BUILTIN,HIGH);
  for (i=(npts>>1);i!=0;i>>=1) bits++;
  fft_init(bits);
  Serial.println("Finished with Setup");
}
//.................................................................loop
void loop() { 
  server2.handleClient();  // handle webserver on port 80
  webSocket.loop();
  if (sample_buffer_ready) {
    digitalWrite(LED_BUILTIN,LOW);
    sample_buffer_ready=0;
    int points_to_send=npts;
    if (spectrum) {
      double qq=0.0;
      for (int i=0;i<npts;i++) { // find average
        qq+=sample_buffer[i]; 
      }
      qq/=npts;
      for (int k=0;k<npts;k++) {
        wave[k].re=sample_buffer[k]-qq;
        wave[k].im=0;          
      }
      yield();
      fft(wave,bits); yield();
      reorder(wave,bits); yield();
      qq=2.0*vscale/npts;
      for (int i=0;i<npts/2;i++) {
        sample_buffer[i]=hypot(wave[i].re,wave[i].im)*qq;
      }
      points_to_send=npts/2;
    } 
    yield();
    doc.to<JsonObject>();  // clear doc
    for (int j=0;j<points_to_send;j++){doc["WF1"][j]=sample_buffer[j];}
    serializeJson(doc,out);
    webSocket.sendTXT(websock_num,out,strlen(out));    
    SampleSlow.attach_ms(time_to_wait,sampleslow_action,websock_num);   
//  Serial.println(json.measureLength()); 
    digitalWrite(LED_BUILTIN,HIGH); 
//  float volt=(analogRead(A0)-4.0)*3.3/(1004.0-4.0);
//  Serial.print("A0 = "); Serial.println(volt);
  }
}
