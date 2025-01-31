#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
int mergeball_kernel(char *fname1, char *fname2, char *fname_out);

int main(void)
{
  char *fname1 = "b31.xyz";
  char *fname2 = "b41.xyz";
  char *fname_out = "mergeball.xyz";

  mergeball_kernel(fname1,fname2,fname_out);

  return 0;
}
*/

int mergeball_kernel(char *fname1, char *fname2, char *fname_out)
{
  const char *my_name = "mergeball_kernel";
  FILE *fd;
  char *str;
  size_t slen;
  char **s1;
  char **s2;
  int i;
  int j;
  int n1;
  int n2;
  int nn2;
  int *uniq;

  // debug
  /*
  printf("fname1 = %s\nfname2 = %s\nfname_out = %s\n",fname1,fname2,fname_out);
  return 33;
  */

  // read files
  str = NULL;

  fd = fopen(fname1,"r");
  if (fd == NULL) {
    printf("%s: cannot open %s\n",my_name,fname1);
    return 1;
  }
  getline(&str,&slen,fd);
  sscanf(str,"%d",&n1);
  free(str);
  str = NULL;
  getline(&str,&slen,fd);
  free(str);
  str = NULL;
  s1 = (char **) malloc(n1 * sizeof(char *));
  for (i = 0; i < n1; i++) {
    getline(&str,&slen,fd);
    s1[i] = str;
    str = NULL;
  }
  fclose(fd);

  fd = fopen(fname2,"r");
  if (fd == NULL) {
    printf("%s: cannot open %s\n",my_name,fname2);
    return 1;
  }
  getline(&str,&slen,fd);
  sscanf(str,"%d",&n2);
  free(str);
  str = NULL;
  getline(&str,&slen,fd);
  free(str);
  str = NULL;
  s2 = (char **) malloc(n2 * sizeof(char *));
  for (i = 0; i < n2; i++) {
    getline(&str,&slen,fd);
    s2[i] = str;
    str = NULL;
  }
  fclose(fd);

  // allocate uniq
  uniq = (int *) malloc(n2 * sizeof(int));
  for (i = 0; i < n2; i++) uniq[i] = 1;

  // compare strings
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      if (strcmp(s1[i],s2[j]) == 0) {
        uniq[j] = 0;
        break;
      }
    }
  }

  // write output
  fd = fopen(fname_out,"w");
  if (fd == NULL) {
    printf("%s: cannot open %s\n",my_name,fname_out);
    return 1;
  }

  nn2 = 0;
  for (i = 0; i < n2; i++) nn2 += uniq[i];

  fprintf(fd,"%d\n\n",n1+nn2);
  for (i = 0; i < n1; i++) {
    fprintf(fd,"%s",s1[i]);
  }
  for (i = 0; i < n2; i++) {
    if (uniq[i]) fprintf(fd,"%s",s2[i]);
  }

  fclose(fd);

  // free memory
  for (i = 0; i < n1; i++) free(s1[i]);
  free(s1);
  for (i = 0; i < n2; i++) free(s2[i]);
  free(s2);
  free(uniq);

  return 0;
}
