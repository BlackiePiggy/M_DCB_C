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
    //P1，P2，L1，L2均为二维double数组
    double **P1;
    double **P2;
    double **L1;
    double **L2;
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
void writeArrayToFileWithLabels(FILE *file, double **array, const char* dataType);


// -----------------------------Main--------------------------------------
int main() {
    printf("elevation angle threshold(unit:degree)(10 degree is recommended): %d\n", lim);
    printf("the order of spheric harmonic function (4 order is recommended): %d\n", order);
    printf("MDCB(multi-stations) starts running!\n");

    // Step one: read rinex files
    SitesInfo sitesInfo;
    Obs obs;
    read_rinex(r_ipath, r_opath, &sitesInfo, &obs);

    printf("Hello, world!\n");

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
    //为obs结构体内的数组分配内存，每个变量均需要一个(2880,32)尺寸的二维double数组
    // 分配内存空间给 P1
    pobs->P1 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        pobs->P1[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            pobs->P1[i][j] = 0.0; // 初始化为零
        }
    }
    // 分配内存空间给 P2
    pobs->P2 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        pobs->P2[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            pobs->P2[i][j] = 0.0; // 初始化为零
        }
    }

    // 分配内存空间给 L1
    pobs->L1 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        pobs->L1[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            pobs->L1[i][j] = 0.0; // 初始化为零
        }
    }

    // 分配内存空间给 L2
    pobs->L2 = (double **)malloc(2880 * sizeof(double *));
    for (int i = 0; i < 2880; i++) {
        pobs->L2[i] = (double *)malloc(32 * sizeof(double));
        for (int j = 0; j < 32; j++) {
            pobs->L2[i][j] = 0.0; // 初始化为零
        }
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
    int obst_n; char **obst;
    double h; double m; double s; int ep = 0;
    int *loc;

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

        //给pobs结构体内的数组赋初值为0
        for (int j = 0; j < 2880; j++) {
            for (int k = 0; k < 32; k++) {
                pobs->P1[j][k] = 0.0;
                pobs->P2[j][k] = 0.0;
                pobs->L1[j][k] = 0.0;
                pobs->L2[j][k] = 0.0;
            }
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
                while (1){
                    //读取下一行
                    fgets(line, sizeof(line), file);
                    // 检查是否到达文件结尾
                    if (feof(file)) {
                        // 到达文件结尾，跳出循环
                        break;
                    }
                    //如果line长度大于32个个字符，并且line的第33个字符是"G"或"R"或空格，则执行下面的语句
                    if (strlen(line) > 32 && (line[32] == 'G' || line[32] == 'R')){
                        //将line的第11到第12个字符转换成double类型，再赋值给h
                        h = 0; m = 0; s = 0;
                        if (sscanf(line + 10, "%2lf", &h) == 1){
                        }
                        if (sscanf(line + 13, "%2lf", &m) == 1){
                        }
                        if (sscanf(line + 16, "%10lf", &s) == 1){
                        }
                        ep = h*120 + m*2 + s/30 + 1;

//                        if (ep == 17){
//                            printf("ep: %d\n", ep);
//                        }
                    }else{
                        continue;
                    }
                    //---get satellite number----

                    //debug专用，如果line的值等于下面的字符串，则打印这行内容
                    int nsat; int *sv_G;
                    //取出line的第31到32个字符，转换成int类型，再赋值给nsat
                    if (sscanf(line + 30, "%2d", &nsat) == 1){
                    }
                    //给sv_G申请nsat个int类型的内存，并将数组中每个值赋值为0
                    sv_G = (int *) malloc(nsat * sizeof(int));
                    for (int j = 0; j < nsat; j++){
                        sv_G[j] = 0;
                    }
                    if (nsat>12){
                        for(int j = 0; j < 12; j++){
                            if (line[29+3*(j+1)] == 'G'){
                                //取出line的第31+3*j到32+3*j个字符，转换成int类型，再赋值给sv_G[j]
                                if (sscanf(line + 30 + 3*(j+1), "%2d", &sv_G[j]) == 1){
                                }
                            }
                        }
                        //读取下一行
                        fgets(line, sizeof(line), file);
                        if (nsat<25){
                            for (int j = 0; j < nsat-12; j++){
                                if (line[29+3*(j+1)] == 'G' || line[32] == ' '){
                                    if (sscanf(line + 30 + 3*(j+1), "%2d", &sv_G[j+12]) == 1){
                                    }
                                }
                            }
                        }
                    }else{
                        for(int j = 0; j < nsat; j++){
                            if (line[29+3*(j+1)] == 'G' || line[32] == ' '){
                                //取出line的第32+3*j到32+3*j+2个字符，转换成int类型，再赋值给sv_G[j]
                                if (sscanf(line + 30 + 3*(j+1), "%2d", &sv_G[j]) == 1){
                                }
                            }
                        }
                    }

                    //---get the observations----
                    for (int j = 0; j < nsat; j++){
                        fgets(line, sizeof(line), file);

                        double *obs_temp;
                        obs_temp = (double *) malloc(obst_n * sizeof(double));
                        //给obs_temp都赋初值为0
                        for (int k = 0; k < obst_n; k++){
                            obs_temp[k] = 0;
                        }

                        if (obst_n>5){
                            if (sv_G[j] == 0){
                                fgets(line, sizeof(line), file);
                                continue;
                            }
                            for (int k = 0; k < 5; k++){//思朗科技
                                //判断line的长度，如果大于16*(k+1)-3，则执行下面的语句
                                if (strlen(line) > 16*(k+1)-3){
                                    //strcmp(line(16*j-15:16*j-2),'              '),如果line的第15*(k+1)-15到第15*(k+1)-2个字符是空格，则执行下面的语句
                                    if (strncmp(line + 16*k, "              ", 14) == 0){
                                        continue;
                                    }
                                    //取出line的第15*(j+1)-15至第15*(j+1)-2个字符，转换成double类型，再赋值给obs_temp[j]
                                    if (sscanf(line + 16*k, "%14lf", &obs_temp[k]) == 1){
//                                        //如果obs_temp[k]的前8位数等于25789317，则打印这行内容
//                                        if (obs_temp[k] == 25789317.7094){
//                                            printf("ep: %d\n", ep);
//                                        }
                                    }
                                }
                            }
                            fgets(line, sizeof(line), file);
                            for (int k = 0; k < obst_n-5; k++) {
                                // 判断line的长度，如果大于16*(k+1)-3，则执行下面的语句
                                if (strlen(line) > 16*(k+1)-3) {
                                    // strcmp(line + 16*k, "              "), 如果line的第16*k至第16*k+14个字符是空格，则执行下面的语句
                                    if (strcmp(line + 16*k, "              ") == 0) {
                                        continue;
                                    }
                                    // 取出line的第16*k至第16*k+14个字符，转换成double类型，再赋值给obs_temp[k+5]
                                    // 注意：C语言中数组索引从0开始，因此对应MATLAB中的j+5需要在C语言中调整为k+5
                                    if (sscanf(line + 16*k, "%14lf", &obs_temp[k+5]) == 1) {
                                        // 成功转换字符串为double并赋值
                                    }
                                }
                            }
                        }else{
                            if (sv_G[j] == 0){
                                fgets(line, sizeof(line), file);
                                continue;
                            }
                            for (int k = 0; k < obst_n; k++){
                                //判断line的长度，如果大于16*(k+1)-3，则执行下面的语句
                                if (strlen(line) > 16*(k+1)-3){
                                    //strcmp(line(16*j-15:16*j-2),'              '),如果line的第15*(k+1)-15到第15*(k+1)-2个字符是空格，则执行下面的语句
                                    if (strncmp(line + 16*k, "              ",14) == 0){
                                        continue;
                                    }
                                    //取出line的第15*(j+1)-15至第15*(j+1)-2个字符，转换成double类型，再赋值给obs_temp[j]
                                    if (sscanf(line + 16*k, "%14lf", &obs_temp[k]) == 1){
                                    }
                                }
                            }
                        }
                        pobs->L1[ep-1][sv_G[j]-1]=obs_temp[loc[0]];
                        pobs->L2[ep-1][sv_G[j]-1]=obs_temp[loc[1]];
                        pobs->P1[ep-1][sv_G[j]-1]=obs_temp[loc[2]];
                        pobs->P2[ep-1][sv_G[j]-1]=obs_temp[loc[3]];
                    }
                }
            }
        }
        // 关闭文件
        fclose(file);

        char sav_filename[256]; // 假设文件名长度不超过256个字符

        sprintf(sav_filename, "D:\\projects\\M_DCB\\RINEX_output_files\\observation_%s.csv", psitesInfo->name[i]);

        // Open file for writing
        FILE *outfile = fopen(sav_filename, "w");
        if (outfile == NULL) {
            fprintf(stderr, "Error opening file for writing\n");
        }

        // Write each array to file with labels
        writeArrayToFileWithLabels(outfile, pobs->P1, "P1");
        writeArrayToFileWithLabels(outfile, pobs->P2, "P2");
        writeArrayToFileWithLabels(outfile, pobs->L1, "L1");
        writeArrayToFileWithLabels(outfile, pobs->L2, "L2");

        // Close file
        fclose(outfile);

        printf("Data written to file successfully.\n");

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

// Function to write array to file with labels
void writeArrayToFileWithLabels(FILE *file, double **array, const char* dataType) {
    // 写入列标签
    fprintf(file, "%s_Epoch", dataType); // 假设第一列为时间戳
    for (int j = 1; j <= 32; j++) {
        fprintf(file, ",%s_sv_%d", dataType, j); // 生成并写入标签
    }
    fprintf(file, "\n"); // 结束标签行

    // 写入数据
    for (int i = 0; i < 2880; i++) {
        fprintf(file, "%d", i+1); // 假设第一列为时间戳，此处简化处理
        for (int j = 0; j < 32; j++) {
            fprintf(file, ",%.4f", array[i][j]);
        }
        fprintf(file, "\n");
    }
}
