#include <stdio.h>
#include <stdlib.h>

void draw(int *** image, int height, int width){
    int x,y; 
    int red,green,blue; 
    FILE *image_file = fopen("image.ppm","w");
    fprintf(image_file, "P6\n"); 
    fprintf(image_file, "%d %d\n",width,height); 
    fprintf(image_file, "255\n"); 
    /* Do this for the whole picture now */ 
    for (y = 0 ; y < height ; y++) { 
	for (x = 0;x < width;x++) {
	    red = image[y][x][0];
	    green = image[y][x][1];
	    blue = image[y][x][2];
	    fputc((char)red,image_file); 
	    fputc((char)green,image_file); 
	    fputc((char)blue,image_file); 
	}
	//		fprintf(image_file, "\n"); 
    }
}


