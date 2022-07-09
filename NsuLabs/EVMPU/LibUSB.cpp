#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <time.h>
#include <stdio.h>
 
int main(int argc,char *argv[])
{ CvCapture *capture = cvCreateCameraCapture(0);
    if (!capture) return 0;
    struct timespec start,end,start1,end1,start2,start3,end2,end3;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    unsigned long long cnt = 0;
    while(1) {
        clock_gettime(CLOCK_MONOTONIC_RAW, &start3);
        IplImage *frame = cvQueryFrame(capture);
        clock_gettime(CLOCK_MONOTONIC_RAW, &end3);
        if(!frame) break;
        clock_gettime(CLOCK_MONOTONIC_RAW, &start1);
        for (int y=0; y < frame->height; y++) {
            uchar *ptr = (uchar*)(frame->imageData + y*frame->widthStep);
            int red,green,blue;
            for (int x=0; x < frame->width; x++) {
                red = ptr[3*x+2];
                green = ptr[3*x+1];
                blue = ptr[3*x];
                ptr[3*x+1] = red; // Green
                ptr[3*x+2] = blue;
                ptr[3*x] = green;
                }
            }
        cvFlip(frame,frame, 1);
        clock_gettime(CLOCK_MONOTONIC_RAW, &end1);
        clock_gettime(CLOCK_MONOTONIC_RAW, &start2);
        cvShowImage("test", frame);
        clock_gettime(CLOCK_MONOTONIC_RAW, &end2);
        char c = cvWaitKey(33);
        cnt+=1;
        if (cnt == 18446744073709551614) break;
        if(c == 27) break;
    }
    cvReleaseCapture(&capture);
    cvDestroyWindow("test");
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Время на преобразование (на 1 кадр) %lf\n",end1.tv_sec-start1.tv_sec + 0.000000001*(end1.tv_nsec-start1.tv_nsec));
    printf("Время на вывод изображения (1 кадра) %lf\n",end2.tv_sec-start2.tv_sec + 0.000000001*(end2.tv_nsec-start2.tv_nsec));
    printf("fps: %f\n",cnt/(end.tv_sec-start.tv_sec + 0.000000001*(end.tv_nsec-start.tv_nsec)));
}
