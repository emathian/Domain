/* strncmp example */
#include <stdio.h>
#include <string.h>

#define LINE_LENGTH     81  


char* substring(const char* str, size_t begin, size_t len);
void test( int* x); 

int main ()
{
  char str[][14] = { "ATOM  R2D2" , "C3PO" , "ATOM  BR2A6" };
  int n;
  puts ("Looking for R2 astromech droids...");

	 char    line[LINE_LENGTH]= "ATOM    182  CB  GLN A  25      17.779  24.176  85.217";
  size_t      begin  = 2; 
  size_t      end    = 4; 

  char*       substr2 = substring(line, begin, end); 
  printf ("print Substr22222  %s\n", substr2);
  //char * str;
  //str = str2.substr(2,4);
  printf ("print %lu\n", strlen(str[0]));
  char*       substr = substring(str[0], begin, end); 
  printf ("print Substr  %s\n", substr);


  if (strncmp (substr,"OM",1) == 0)
    {
      printf ("found sub %s\n",substr);
    }


  for (n=0 ; n<3 ; n++)
    if (strncmp (str[n],"ATOM  xR2xx",8) == 0)
    {
      printf ("found %s\n",str[n]);
    }

    int X = 0;
    test(&X);
    printf("%d\n", X);


    int *A;
    &A = 3;
    printf("%d\n", &A);

  return 0;
}

char* substring(const char* str, size_t begin, size_t len) 
{ 
  if (str == 0 || strlen(str) == 0 || strlen(str) < begin || strlen(str) < (begin+len)) 
    return 0; 

  return strndup(str + begin, len); 
} 


void test( int* x) 
{ 
  *x = 5;
} 