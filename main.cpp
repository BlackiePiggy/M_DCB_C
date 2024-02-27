#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>


// 假设 SitesInfo 结构体
typedef struct {
    char** name;
    int* doy;
    double** coor;
} SitesInfo;

typedef struct {
    double* P1;
    double* P2;
    double* L1;
    double* L2;
} Obs;

// -----------------------------Setting--------------------------------------
const char* r_ipath = "D:\\projects\\M_DCB\\RINEX_files";// rinex
const char* s_ipath = "D:\\projects\\M_DCB\\SP3_files";// sp3
const char* i_ipath = "D:\\projects\\M_DCB\\IONEX_files";// ionex
const char* r_opath = "M_OBS";
int lim = 10;// el
int order = 4;// order

// -----------------------------Function--------------------------------------
int countFilesInDirectory(const char *folderPath);
//void free_usage(SitesInfo sitesInfo, int len);
void read_rinex(const char* r_ipath, const char* r_opath, SitesInfo *psitesInfo, Obs *pobs);

// -----------------------------Main--------------------------------------
int main() {
    printf("elevation angle threshold(unit:degree)(10 degree is recommended): %d\n", lim);
    printf("the order of spheric harmonic function (4 order is recommended): %d\n", order);
    printf("MDCB(multi-stations) starts running!\n");

    // Step one: read rinex files
    SitesInfo sitesInfo;
    Obs obs;
    read_rinex(r_ipath, r_opath, &sitesInfo, &obs);


    // 释放 sitesInfo 结构体内存
    //free_usage(sitesInfo, len);

    return 0;
}

void read_rinex(const char* r_ipath, const char* r_opath, SitesInfo *psitesInfo, Obs *pobs) {
    DIR *dir;
    FILE *file;
    struct dirent *entry;
    int counter;
    int r_file_num;

    // 打开文件夹
    dir = opendir(r_ipath);
    if (dir == NULL) {
        perror("无法打开文件夹");
        exit(EXIT_FAILURE);
    }

    printf("Starting Reading Rinex Files!\n");

    // 计算文件数
    r_file_num = countFilesInDirectory(r_ipath);

    // 为 SitesInfo 结构体内的字符串数组 name 分配内存
    psitesInfo->name = (char **) malloc(r_file_num * sizeof(char *));
    for (int i = 0; i < r_file_num; i++) {
        psitesInfo->name[i] = (char *) malloc(256 * sizeof(char)); // 假设文件名最大长度为256
    }
    psitesInfo->doy = (int *) malloc(r_file_num * sizeof(int));
    psitesInfo->coor = (double **) malloc(r_file_num * sizeof(double *));
    for (int i = 0; i < r_file_num; i++) {
        psitesInfo->coor[i] = (double *) malloc(256 * sizeof(double)); // 假设文件名最大长度为256
    }


    int index = 0;
    int i = 0;
    char doy[4];
    char sub_str[3];
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[index] != '.') {
            strcpy(psitesInfo->name[i], entry->d_name);
            strncpy(sub_str, &entry->d_name[9], 2);
            strncpy(doy, &entry->d_name[4], 3);
            sub_str[2] = '\0';
            doy[3] = '\0';  //这一步必须加入字符串标志，否则会出问题
            //printf("sub_str:%d\n", atoi(sub_str));
            //printf("doy:%d\n", atoi(doy));
            psitesInfo->doy[i] = 1000 * atoi(sub_str) + atoi(doy);
            i++;
        }
    }

    char line[256]; // 假设一行最多包含256个字符
    char filename[256];

    for (int i = 0; i < r_file_num; i++) {
        printf("Reading File No.%d :%s\n", i + 1, psitesInfo->name[i]);
        // 打开文件
        strcpy(filename, r_ipath);
        strcat(filename, "\\");
        strcat(filename, psitesInfo->name[i]);
        file = fopen(filename, "r");
        if (file == NULL) {
            perror("Can't Open File");
        }

        while (fgets(line, sizeof(line), file)) {
            // -----get receivers coordinates-----------
            if (strstr(line, "APPROX POSITION XYZ") != NULL) {
                //printf("%s", line);
                double x, y, z;
                // 将行解析成三个双精度浮点数
                if (sscanf(line, "%lf %lf %lf", &x, &y, &z) == 3) {
                    psitesInfo->coor[i][0] = x;
                    psitesInfo->coor[i][1] = y;
                    psitesInfo->coor[i][2] = z;
                }
            }
            //-----get the GPS observables' types------
            if (strstr(line, "# / TYPES OF OBSERV") != NULL){
                //printf("%s", line);
                int obst_n; char **obst;
                if(sscanf(line, "%d", &obst_n) == 1){
                    if (obst_n > 9){
                        printf("obst_n > 9!!!!!");
                    }
                    else{
                        obst = (char **) malloc(obst_n * sizeof(char *));
                        for (int j = 0; j < obst_n; j++) {
                            obst[j] = (char *) malloc(4 * sizeof(char)); // 假设每个观测量最大长度为4
                        }
                        int offset = 10;
                        for (int j = 0; j < obst_n; j++){
                            if (sscanf(line + offset, "%2s", obst[j]) == 1) {
                                // 成功解析一个观测名称，更新偏移量，准备解析下一个观测名称
                                offset += 6; // 每个观测名称占据4个字符，包括空格
                            } else {
                                // 解析失败，处理异常情况
                                printf("Failed to parse observation type at index %d\n", j);
                                break;
                            }
                        }
                        int *loc;
                        loc = (int *) malloc(obst_n * sizeof(int));
                        for (int j = 0; j < obst_n; j++){
                            //比较obst[j]和"L1"字符串是否相同，如果相同，则给loc[0]赋值为j
                            if (strcmp(obst[j], "L1") == 0){
                                loc[0] = j;
                            }else if (strcmp(obst[j], "L2") == 0) {
                                loc[1] = j;
                            }else if (strcmp(obst[j], "P1") == 0) {
                                loc[2] = j;
                            }else if (strcmp(obst[j], "P2") == 0) {
                                loc[3] = j;
                            }
                        }
                        if (loc[2] == 0 || loc[3] == 0){
                            printf("P1 and P2 not found!\n");
                        }
                    }
                }
            }
            //----start get GPS observables------------
            if (strstr(line, "END OF HEADER") != NULL) {
                while (fgets(line, sizeof(line), file)){
                    //如果line长度大于32个个字符，并且line的第33个字符是"G"或"R"或空格，则执行下面的语句
                    if (strlen(line) > 32 && (line[32] == 'G' || line[32] == 'R' || line[32] == ' ')){
                        //将line的第11到第12个字符转换成double类型，再赋值给h
                        double h; double m; double s; double ep;
                        if (sscanf(line + 11, "%1lf", &h) == 1){
                        }
                        if (sscanf(line + 14, "%1lf", &m) == 1){
                        }
                        if (sscanf(line + 17, "%10lf", &s) == 1){
                        }
                        ep = h*120 + m*2 + s/30 + 1;
                    }
                }
            }
        }
        // 关闭文件
        fclose(file);
    }
    printf("Step one: completing !\n");
}

// 文件夹内文件数
int countFilesInDirectory(const char *folderPath) {
    DIR *dir;
    struct dirent *entry;
    int counter;

    // 打开文件夹
    dir = opendir(folderPath);

    if (dir == NULL) {
        perror("无法打开文件夹");
        exit(EXIT_FAILURE);
    }

    // 读取文件夹中的文件
    while ((entry = readdir(dir)) != NULL) {
        counter ++;
//        printf("No.%d :%s\n", counter,entry->d_name); //Print file names
    }
    counter = counter - 2;

    // 关闭文件夹
    closedir(dir);

    return counter;
}

//void free_usage(SitesInfo sitesInfo, int len) {
//    // 释放 sitesInfo 结构体内存
//    for (int i = 0; i < len; i++) {
//        free(sitesInfo.name[i]);
//    }
//    free(sitesInfo.name);
//    free(sitesInfo.doy);
//    for (int i = 0; i < len; i++) {
//        free(sitesInfo.coor[i]);
//    }
//    free(sitesInfo.coor);
//}