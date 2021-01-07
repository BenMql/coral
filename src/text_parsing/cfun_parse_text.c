#include <stdio.h>
#include <string.h>

int parse_text(int *numberOfLines, int *phase, int *lengths, char (*arr_lines)[1024] )
{
    FILE* file = fopen("coral.equations", "r"); /* should check the result */
    char line[1024];
    char *trimEOLchar;

    int current_line=0;

    while (fgets(line, sizeof(line), file)) {
	current_line+=1;
        trimEOLchar = strtok(line, "\n");
	if (*phase > 1) {
	lengths[current_line-1] = strlen(line);
        strcpy(arr_lines[current_line-1], line);
	}
    }

    *numberOfLines=current_line;
    fclose(file);

    return 0;
}
int parse_output(int *numberOfLines, int *phase, int *lengths, char (*arr_lines)[1024] )
{
    FILE* file = fopen("coral.usrOutput", "r"); /* should check the result */
    char line[1024];
    char *trimEOLchar;

    int current_line=0;

    while (fgets(line, sizeof(line), file)) {
	current_line+=1;
        trimEOLchar = strtok(line, "\n");
	if (*phase > 1) {
	lengths[current_line-1] = strlen(line);
        strcpy(arr_lines[current_line-1], line);
	}
    }

    *numberOfLines=current_line;
    fclose(file);

    return 0;
}

int parse_timeseries(int *numberOfLines, int *phase, int *lengths, char (*arr_lines)[1024] )
{
    FILE* file = fopen("coral.timeseries", "r"); /* should check the result */
    char line[1024];
    char *trimEOLchar;

    int current_line=0;

    while (fgets(line, sizeof(line), file)) {
	current_line+=1;
        trimEOLchar = strtok(line, "\n");
	if (*phase > 1) {
	lengths[current_line-1] = strlen(line);
        strcpy(arr_lines[current_line-1], line);
	}
    }

    *numberOfLines=current_line;
    fclose(file);

    return 0;
}
