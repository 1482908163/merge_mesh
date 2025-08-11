#include <iostream>
#include <climits>
#include "mpi.h"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Shape.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "3DNgmesher.h"
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

namespace nglib {
#include <nglib.h>
}
#include "createElmerOutput.h"

void print_help() {
    cerr << "参数输入帮助" << endl <<
         "-o : 输出文件目录 默认为./output" << endl <<
         "-i : 输入文件目录 默认为/work/home/moussa/wholewall3.stp" << endl <<
         "-l --levels : 细化等级 默认为0" << endl <<
         "-r --refine : 细化次数 默认为0" << endl <<
         "-maxh : 网格最大值 默认为1000.0" << endl <<
         "-minh : 网格最小值 默认为10.0" << endl <<
         "-v : 保存细化文件" << endl <<
         "-adj : 通信" << endl <<
         "-h --help : 参数输入帮助" << endl;
}

int main(int argc, char **argv) {

    using namespace nglib;

    int id;
    int p = 1;
    MPI_Init(nullptr, nullptr);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // if(id == 0) cout << "MPI:" << p << endl;


    Ng_Init();

    Ng_Mesh *occ_mesh;

//Parameters
    string OUTPUT_PATH = "./output/";
    string INPUT_PATH = "/work/home/moussa/wholewall3.stp";
    bool save_vol = false;
    bool isComputeAdj = false;
    int numlevels = 0;
    int numrefine = 0;
    double maxh = 1000.0,minh = 10.0;
//

    if(argc <= 1) {
        print_help();
        MPI_Finalize();
        return 1;
    }

    for(int i = 1; i < argc; i++) {
        if(!strcmp(argv[i],"-o")) {
            if(argv[i+1] != NULL) OUTPUT_PATH = argv[i+1];
            else {
                print_help();
                MPI_Finalize();
                return 1;
            }
        }
        else if(!strcmp(argv[i],"-i")) {
            if(argv[i+1] != NULL) INPUT_PATH = argv[i+1];
            else {
                print_help();
                MPI_Finalize();
                return 1;
            }
        }
        else if(!strcmp(argv[i],"-l") || !strcmp(argv[i],"--levels")) {
            if(argv[i+1] != NULL) numlevels = atoi(argv[i+1]);
            else {
                print_help();
                MPI_Finalize();
                return 1;
            }
        }
        else if(!strcmp(argv[i],"-r") || !strcmp(argv[i],"--refine")) {
            if(argv[i+1] != NULL) {
                numrefine = atoi(argv[i+1]);
            }
            else {
                print_help();
                MPI_Finalize();
                return 1;
            }
        }
        else if(!strcmp(argv[i],"--maxh")) {
            if(argv[i+1] != NULL) maxh = atof(argv[i+1]);
            else {
                print_help();
                MPI_Finalize();
                return 1;
            }
        }
        else if(!strcmp(argv[i],"--minh")) {
            if(argv[i+1] != NULL) minh = atof(argv[i+1]);
            else {
                print_help();
                MPI_Finalize();
                return 1;
            }
        }
        else if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
            print_help();
            MPI_Finalize();
            return 1;
        }
        else if(!strcmp(argv[i],"-v")) {
            save_vol = true;
        }
        else if(!strcmp(argv[i],"-adj")) {
            isComputeAdj = true;

        }
    }

    string savepvname = OUTPUT_PATH + "test_occ/test_occ.vol";

    if(id == 0) {

        cerr << "========参数========" << endl <<
             "使用核数 : " << p << endl <<
             "输出文件目录 : " << OUTPUT_PATH << endl <<
             "输入文件目录 : " << INPUT_PATH << endl <<
             "细化等级 : " << numlevels << endl <<
             "细化次数 : " << numrefine << endl <<
             "网格最大值:" << maxh << endl <<
             "网格最小值:" << minh << endl;


        // int ret = mkdir(OUTPUT_PATH.c_str(), 0777);
        mkdir(OUTPUT_PATH.c_str(), 0777);

        string meshQuality_path = OUTPUT_PATH + string("meshQuality/");
        // int meshQuality_ret = mkdir(meshQuality_path.c_str(), 0777);
        mkdir(meshQuality_path.c_str(), 0777);

        string refinedSurfmesh_path = OUTPUT_PATH + string("refinedSurfmesh/");
        // int refinedSurfmesh_ret = mkdir(refinedSurfmesh_path.c_str(), 0777);
        mkdir(refinedSurfmesh_path.c_str(), 0777);

        string test_occ_path = OUTPUT_PATH + string("test_occ/");
        // int test_occ_ret = mkdir(test_occ_path.c_str(), 0777);
        mkdir(test_occ_path.c_str(), 0777);

        string testout_path = OUTPUT_PATH + string("testout/");
        // int testout_ret = mkdir(testout_path.c_str(), 0777);
        mkdir(testout_path.c_str(), 0777);

        string volfined_path = OUTPUT_PATH + string("volfined/");
        // int volfined_ret = mkdir(volfined_path.c_str(), 0777);
        mkdir(volfined_path.c_str(), 0777);

        string volwithadj_path = OUTPUT_PATH + string("volwithadj/");
        // int volwithadj_ret = mkdir(volwithadj_path.c_str(), 0777);
        mkdir(volwithadj_path.c_str(), 0777);
    }

    // Define pointer to OCC Geometry
    Ng_OCC_Geometry *occ_geom;

    // Ng_Mesh *occ_mesh;

    Ng_Meshing_Parameters mp;

    TopTools_IndexedMapOfShape FMap;

    // Ng_OCC_TopTools_IndexedMapOfShape *occ_fmap = (Ng_OCC_TopTools_IndexedMapOfShape*)&FMap;

    // Result of Netgen Operations
    Ng_Result ng_res;

    // Initialise the Netgen Core library
    // Ng_Init();

    // Read in the OCC File
    MPI_Barrier(MPI_COMM_WORLD);
    double startTime = MPI_Wtime();

    string STEP_PATH = INPUT_PATH;
    occ_geom = Ng_OCC_Load_STEP(STEP_PATH.c_str());
    if (!occ_geom)
    {
        cout << "Error reading in STEP File: " << STEP_PATH << endl;
        MPI_Finalize();
        return 1;
    }
    if(id == 0)
        cout << "Successfully loaded STEP File: " << STEP_PATH << endl;


    occ_mesh = Ng_NewMesh();

    mp.uselocalh = 1;
    mp.elementsperedge = 2.0;
    mp.elementspercurve = 2.0;
    mp.maxh = maxh;
    mp.minh = minh;
    mp.grading = 0.3;
    mp.closeedgeenable = 0;
    mp.closeedgefact = 1.0;
    mp.optsurfmeshenable = 0;


    if(id == 0) {
        cout << "Setting Local Mesh size....." << endl;
        cout << "OCC Mesh Pointer before call = " << occ_mesh << endl;
    }
    Ng_OCC_SetLocalMeshSize(occ_geom, occ_mesh, &mp);
    if(id == 0) {
        cout << "Local Mesh size successfully set....." << endl;
        cout << "OCC Mesh Pointer after call = " << occ_mesh << endl;
        cout << "Creating Edge Mesh....." << endl;
    }

    ng_res = Ng_OCC_GenerateEdgeMesh(occ_geom, occ_mesh, &mp);
    if (ng_res != NG_OK)
    {
        Ng_DeleteMesh(occ_mesh);
        cout << "Error creating Edge Mesh.... Aborting!!" << endl;
        MPI_Finalize();
        return 1;
    }
    else
    {
        if(id == 0) {
            cout << "Edge Mesh successfully created....." << endl;
            cout << "Number of points = " << Ng_GetNP(occ_mesh) << endl;
        }
    }

    id == 0 ? cout << "Creating Surface Mesh....." << endl: cout << "";

    ng_res = Ng_OCC_GenerateSurfaceMesh(occ_geom, occ_mesh, &mp);
    if (ng_res != NG_OK)
    {
        Ng_DeleteMesh(occ_mesh);
        cout << "Error creating Surface Mesh..... Aborting!!" << endl;
        MPI_Finalize();
        return 1;
    }
    else
    {
        if(id == 0) {
            cout << "Surface Mesh successfully created....." << endl;
            cout << "Number of points = " << Ng_GetNP(occ_mesh) << endl;
            cout << "Number of surface elements = " << Ng_GetNSE(occ_mesh) << endl;
        }
    }

    if(id == 0)
        cout << "Creating Volume Mesh....." << endl;

    ng_res = Ng_GenerateVolumeMesh(occ_mesh, &mp);

    if(id == 0) {

        cout << "Volume Mesh successfully created....." << endl;
        cout << "Number of points = " << Ng_GetNP(occ_mesh) << endl;
        cout << "Number of volume elements = " << Ng_GetNE(occ_mesh) << endl;
    }


    // if(id == 0) cout << "Saving Mesh as VOL file....." << endl;
    int faceNum = nglib::Ng_GetNFD((nglib::Ng_Mesh*)occ_mesh);
    int x[4];
    for (int i = 1; i < faceNum; i++) {
        nglib::My_Ng_GetFaceDescriptor((nglib::Ng_Mesh*)occ_mesh, i, x);
        //std::cout <<"+++++++++++++++" << x[0] << "=" << x[1] << "=" << x[2] << "=" << x[3] << std::endl;

    }

    if(id == 0)
        Ng_SaveMesh(occ_mesh, savepvname.c_str());


    if(id == 0) cout << "Generate Coarse Mesh Done..." << endl;
    MPI_Barrier(MPI_COMM_WORLD);


    // double Coarse_endTime = clock();
    double Coarse_endTime = MPI_Wtime();
    double Coarse_Time = (double)(Coarse_endTime - startTime);



    if (p >= 1) {

        double time[5] = {0,0,0,0,0};
        //Set the number of partitions
        int numParts = p;
        FILE *fp;
        FILE *fp_time;
        //set the level of refinement
        // int optvolmeshenable = 0;
        // int optsteps_3d = 0;
        // double gradingp = 0.3;
        // char paramname[20];

        //maxbarycoord is the n secondary of 2
        int maxbarycoord = 1 << (numlevels + numrefine+1);
        //int maxbarycoord = 1 << (numlevels + numrefine);
        map< int, xdMeshFaceInfo > facemap;
        map< int, int > g2lvrtxmap;
        map< Barycentric, int, CompBarycentric > baryc2locvrtxmap;
        map< int, list<int>> adjbarycs;
        list<xdFace> newfaces;
        list<VEindex> VEindexs;
        list<VEindex>::iterator VEi;
        std::map< IntPair, int, IntPairCompare > edgemap;

        std::string str_id = std::to_string(id);

        nglib::Ng_Mesh * submesh = nglib::Ng_NewMesh();
        NewSubmesh(occ_mesh, submesh);

        idx_t *edest = PartitionMesh(occ_mesh, numParts);

        //MPI_Barrier(MPI_COMM_WORLD);
        double currtime0 = MPI_Wtime();
        time[0] = double(currtime0 - Coarse_endTime);

        ExtractPartitionSurfaceMesh(occ_mesh, edest, facemap);
        //MPI_Barrier(MPI_COMM_WORLD);
        double currtime1 = MPI_Wtime();
        time[1] = double(currtime1 - currtime0);


        PartFaceCreate(occ_mesh, id, facemap, maxbarycoord, submesh, g2lvrtxmap, baryc2locvrtxmap, newfaces);
        // savepvname = OUTPUT_PATH + "PartFaceCreate/PartFaceCreate" + str_id + ".vol";
        // if(save_vol) {
        // 	Ng_SaveMesh(submesh, savepvname.c_str());
        // }

        //MPI_Barrier(MPI_COMM_WORLD);
        double currtime2 = MPI_Wtime();
        time[2] = double(currtime2 - currtime1);

        Refine(submesh, numlevels, id, newfaces, baryc2locvrtxmap, edgemap);
        //MPI_Barrier(MPI_COMM_WORLD);
        double currtime3 = MPI_Wtime();
        time[3] = double(currtime3 - currtime2);


        savepvname = OUTPUT_PATH + "refinedSurfmesh/refinedSurfmesh" + str_id + ".vol";
        if(save_vol) {
            Ng_SaveMesh(submesh, savepvname.c_str());
        }




        //The face mesh grid is refined in parallel to each partition.


        int i;

        nglib::Ng_Meshing_Parameters nmp;
        //nmp.maxh = 1e6;
        nmp.fineness = 1;

        double volumeMesh_start = MPI_Wtime();
        nglib::Ng_GenerateVolumeMesh(submesh, &nmp);
        double volumeMesh_end = MPI_Wtime();
        //MPI_Barrier(MPI_COMM_WORLD);
        if(id == 0) printf("meshing done \n");
        double currtime4 = MPI_Wtime();
        time[4] = double(currtime4 - currtime3);

        for (i = 0; i < numrefine; i++) {
            Refineforvol(submesh, id, newfaces, baryc2locvrtxmap, edgemap);
        }

        std::string savepvname = OUTPUT_PATH + "volfined/volfined" + str_id + ".vol";

        if(save_vol) {
            nglib::Ng_SaveMesh(submesh, savepvname.c_str());
        }


        if (isComputeAdj) {
            map<Barycvrtx, list<int>, CompBarycvrtx> barycvrtx2adjprocsmap;
            computeadj(id,facemap,g2lvrtxmap, barycvrtx2adjprocsmap);

            int *VEgid;
            int numNEs = nglib::Ng_GetNE(submesh);
            MYCALLOC(VEgid, int *, (numNEs + 1), sizeof(int));

            printf("start com_barycoords, id: %d\n", id);
            // cout << id << "start com_barycoords" << endl;

            int *newid = com_barycoords(submesh, MPI_COMM_WORLD, barycvrtx2adjprocsmap,
                                        baryc2locvrtxmap, adjbarycs, numParts, VEgid, id);

            
            // int pointdebug = nglib::Ng_GetNP((nglib::Ng_Mesh *)submesh);
            // char *debugpath = new char[512];
            // sprintf(debugpath,"pointdebugpath%d.txt",id);
            // std::string pointdebugpath = OUTPUT_PATH + debugpath;
            // ofstream outpointdebugpath1(pointdebugpath.c_str());
            // for(int i=1;i<=pointdebug;i++){
            //     outpointdebugpath1 << newid[i] << endl;
            // }
            // outpointdebugpath1.close();


            //newid 本地点id--》全局点id
            //VEgid 本地体网格id --》全局体网格id
            //adjbarycs 共享点再哪些处理器上

            // try{
            //     createElmerOutput(submesh,VEgid,newid,adjbarycs,numParts,id,OUTPUT_PATH);
            // }catch(...){
            //     printf("createElmerOutput error id is %d\n", id);
            // }
#if 1
            nglib::Ng_Mesh* mesh = (nglib::Ng_Mesh *)submesh;
            char * boundaryfile1 = new char[512];
            char * elementfile1 = new char[512];
            char * headerfile1 = new char[512];
            char * nodefile1 = new char[512];
            char * sharedfile1 = new char[512];
            char * path1 = new char[512];

            sprintf(path1,"partitioning.%d",numParts);
            sprintf(boundaryfile1, "partitioning.%d/part.%d.boundary", numParts, id+1);
            sprintf(elementfile1, "partitioning.%d/part.%d.elements", numParts, id+1);
            sprintf(headerfile1, "partitioning.%d/part.%d.header", numParts, id+1);
            sprintf(nodefile1, "partitioning.%d/part.%d.nodes", numParts, id+1);
            sprintf(sharedfile1, "partitioning.%d/part.%d.shared", numParts, id+1);

            string path = OUTPUT_PATH + string(path1);
            string boundaryfile = OUTPUT_PATH + string(boundaryfile1);
            string elementfile = OUTPUT_PATH + string(elementfile1);
            string headerfile = OUTPUT_PATH + string(headerfile1);
            string nodefile = OUTPUT_PATH + string(nodefile1);
            string sharedfile = OUTPUT_PATH + string(sharedfile1);

            mkdir(path.c_str(),0777);

            int ne = nglib::Ng_GetNE(mesh); //体网格的数量
            int nse = nglib::Ng_GetNSE(mesh); //面网格的数量
            int np = nglib::Ng_GetNP(mesh); //点的数量

            //输出elements文件
            ofstream outelements(elementfile.c_str());
            int tet[4];
            for(int i=0;i < ne;i++){
                nglib::Ng_GetVolumeElement (mesh, i+1, tet);
                outelements << VEgid[i+1] << " 1 504 " << newid[tet[0]] << " " << newid[tet[1]] << " " << newid[tet[2]] << " " << newid[tet[3]] << endl;
            }
            outelements.close();

            //输出nodes文件
            ofstream outnodes(nodefile.c_str());
            double point[3];
            for(int i=0; i<np;i++){
                nglib::Ng_GetPoint (mesh, i+1, point);
                outnodes << newid[i+1] << " -1 " << point[0] << " " << point[1] << " " << point[2] << endl;
            }
            outnodes.close();

            //输出shared文件
            ofstream outshareds(sharedfile.c_str());
            for(auto it = adjbarycs.begin();it != adjbarycs.end();it++){
                int locid = it->first;
                std::list<int> proceid = it->second;
                int sizeid = proceid.size() + 1;
                std::string sharednode;
                sharednode += std::to_string(sizeid);
                sharednode += " ";
                sharednode += std::to_string(id+1); //elmerID = 核心ID + 1
                sharednode += " ";
                int m = 0;
                for(auto itt = proceid.begin(); itt != proceid.end() ; itt++){
                    sharednode += std::to_string(((*itt)+1));     //elmerID = 核心ID + 1
                    if( m != (sizeid-1) ){
                        sharednode += " ";
                    }
                    m++;
                }
                outshareds << newid[locid] << " " << sharednode << endl;
            }
            outshareds.close();


            //输出boundary文件
            ofstream outboundarys(boundaryfile.c_str());
            //求边界面网格所在的体网格
            Index3 i3;
            int l;
            bool (*fn_pt)(Index3,Index3) = fncomp;
            std::multimap<Index3,int, bool(*)(Index3, Index3)> face2vol(fn_pt);
            std::multimap<Index3,int, bool(*)(Index3, Index3)>::iterator myit;
            for(int i=1; i<=ne;i++){
                nglib::Ng_GetVolumeElement (mesh, i, tet);
                for (int j = 1; j <= 4; j++){
                    l = 0;
                    for (int k = 1; k <= 4; k++)
                    {
                        if (k != j)
                        {
                            i3.x[l] = newid[tet[k-1]];
                            l++;
                        }
                    }
                    i3.Sort();
                    face2vol.insert(pair<Index3,int>(i3,VEgid[i]));
                }
            }

            int *surfpointss = new int[3];
            int surfidx;
            int geoid;
            int number = 0;
            //int nfd = ((Mesh*)mesh)->GetNFD();
            for(int j=0; j< nse; j++){
                nglib::Ng_GetSurfaceElement(mesh, j + 1, surfpointss, surfidx);
                // if((Mesh*)mesh->GetFaceDesriptor(mesh->SurfaceElement(j).GetIndex()).BCProperty()==nfd){
                //     continue;
                // }
                geoid = nglib::GetBoundaryID(mesh,j+1) +1;
                if(nglib::ispatbound(mesh,j+1)){
                    continue;
                }
                i3.x[0] = newid[surfpointss[0]];
                i3.x[1] = newid[surfpointss[1]];
                i3.x[2] = newid[surfpointss[2]];
                i3.Sort();
                myit = face2vol.find(i3);
                if(myit!= face2vol.end()){
                    number++;
                    outboundarys << number << " " << geoid << " " << myit->second << " " << "0" << " 303 " << newid[surfpointss[0]]  << " " << newid[surfpointss[1]] << " " << newid[surfpointss[2]] <<endl;
                }
            }
            outboundarys.close();

            //输出header文件
            ofstream outheader(headerfile.c_str());
            outheader << np << " " << ne << " " << number << endl;
            outheader << 2 << endl;
            outheader << "504 " << ne << endl;
            outheader << "303 " << number << endl;
            if(adjbarycs.size() != 0)
            {
                outheader << adjbarycs.size() << " 0" << endl;
            }
            outheader.close();


#endif

            printf("start com_baryVolumeElements, id: %d\n", id);
            com_baryVolumeElements(submesh, MPI_COMM_WORLD, barycvrtx2adjprocsmap,
                                   baryc2locvrtxmap, adjbarycs, newid, VEgid, VEindexs, numParts, id);
            printf("createElmerOutput, id: %d\n", id);
            

            savepvname = OUTPUT_PATH + "volwithadj/volwithadj" + str_id + ".vol";
            if(save_vol) {
                nglib::Ng_SaveMesh(submesh, savepvname.c_str());

                // string openfoampath = OUTPUT_PATH + "openfoam/part" + str_id;
                // mkdir(openfoampath.c_str(), 0777);
                // const std::filesystem::path  &outfile = openfoampath;
                // nglib::My_WriteOpenFOAMFormat(submesh,outfile);
            }


        }

        /*double endTime = MPI_Wtime();
        double Fine_Time = (double)(endTime - Coarse_endTime);
        double runtime = (double)(endTime - startTime);

        savepvname = OUTPUT_PATH + "volwithadj/volwithadj" + str_id + ".vol";
        if(save_vol) {

        nglib::Ng_SaveMesh(submesh, savepvname.c_str());
        }
        if(id == 0) {
            string savepvname_time = OUTPUT_PATH + "testout/testout_time" + str_id + ".txt";
            fp_time = fopen(savepvname_time.c_str(), "w");
            if (fp_time == NULL) {
                cout << "File " << savepvname << "canot open" << endl;
            }
            else {
            fprintf(fp_time, "Coarse_Time for id:%d is %.2f s\r\n", id, Coarse_Time);
            fprintf(fp_time, "Fine_Time for id:%d is %.2f s\r\n", id, Fine_Time);
            fprintf(fp_time, "runtime for id:%d is %.2f s\r\n", id, runtime);
            for(int i = 0; i < 5; i++) {
                fprintf(fp_time, "part %d time : %.2f \n", i, time[i]);
            }
            }
            fclose(fp_time);
        }

        savepvname = OUTPUT_PATH + "testout/testout_mesh.txt";
        fp = fopen(savepvname.c_str(), "a");
        if (fp == NULL) {
            cout << "File " << savepvname << "canot open" << endl;
        }
        else {
            //fprintf(fp, "the num of points for id:%d is %d\r\n", id, nglib::Ng_GetNP(submesh));
            //fprintf(fp, "the num of Surelemments for id:%d is %d\r\n", id, nglib::Ng_GetNSE(submesh));
            fprintf(fp, "the num of Volelements for id:%d is %d\r\n", id, nglib::Ng_GetNE(submesh));
            fprintf(fp, "the volmesh generate time for if: %d is %f\r\n", id, volumeMesh_end-volumeMesh_start);
        }
        if(id == 0) {
            int Volelements_Sum = nglib::Ng_GetNE(submesh);
            for(int i = 1; i < p; i++) {
                int Volelements_Buf = 0;
                MPI_Recv(&Volelements_Buf, sizeof(Volelements_Buf), MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Volelements_Sum += Volelements_Buf;
            }
            fprintf(fp, "the Sum of Volelements id %d\r\n", Volelements_Sum);

        }
        else {
            int Volelements_Buf = nglib::Ng_GetNE(submesh);
            MPI_Send(&Volelements_Buf, sizeof(Volelements_Buf), MPI_INT, 0, 0, MPI_COMM_WORLD);
        }


        fclose(fp);		/*
        int *VEgids_list, *VEgid_isin_list;
        MYCALLOC(VEgids_list, int *, (VEindexs.size() + 1), sizeof(int));
        MYCALLOC(VEgid_isin_list, int *, (VEindexs.size() + 1), sizeof(int));
        i = 1;
        for (VEi = VEindexs.begin(); VEi != VEindexs.end(); ++VEi) {
            VEgids_list[i] = (*VEi).gid;
            VEgid_isin_list[i] = (*VEi).Isin;
            i++;
        }
        */
        //}
        //else {
        double endTime = MPI_Wtime();
        double Fine_Time = (double)(endTime - Coarse_endTime);
        double runtime = (double)(endTime - startTime);
        savepvname = OUTPUT_PATH + "testout/testout_mesh.txt";
        if(id == 0) {
            string savepvname_time = OUTPUT_PATH + "testout/testout_time" + str_id + ".txt";
            fp_time = fopen(savepvname_time.c_str(), "w");
            if (fp_time == NULL) {
                cout << "File " << savepvname << "canot open" << endl;
            }
            else {
                fprintf(fp_time, "Coarse_Time for id:%d is %.2f s\r\n", id, Coarse_Time);
                fprintf(fp_time, "Fine_Time for id:%d is %.2f s\r\n", id, Fine_Time);
                fprintf(fp_time, "runtime for id:%d is %.2f s\r\n", id, runtime);
                for(int i = 0; i < 5; i++) {
                    fprintf(fp_time, "part %d time : %.2f \n", i, time[i]);
                }
            }
            fclose(fp_time);
        }
        // savepvname = OUTPUT_PATH + "testout/testout.txt";

        fp = fopen(savepvname.c_str(), "a");
        if (fp == NULL) {
            cout << "File " << savepvname << "canot open" << endl;
        }
        else {
            //fprintf(fp, "the num of points for id:%d is %d\r\n", id, nglib::Ng_GetNP(submesh));
            //fprintf(fp, "the num of Surelemments for id:%d is %d\r\n", id, nglib::Ng_GetNSE(submesh));
            fprintf(fp, "the num of Volelements for id:%d is %d\r\n", id, nglib::Ng_GetNE(submesh));
            fprintf(fp, "the volmesh generate time for if: %d is %f\r\n", id, volumeMesh_end-volumeMesh_start);
        }
        if(id == 0) {
            int Volelements_Sum = nglib::Ng_GetNE(submesh);
            for(int i = 1; i < p; i++) {
                int Volelements_Buf = 0;
                MPI_Recv(&Volelements_Buf, sizeof(Volelements_Buf), MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Volelements_Sum += Volelements_Buf;
            }
            fprintf(fp, "the Sum of Volelements id %d\r\n", Volelements_Sum);

        }
        else {
            int Volelements_Buf = nglib::Ng_GetNE(submesh);
            MPI_Send(&Volelements_Buf, sizeof(Volelements_Buf), MPI_INT, 0, 0, MPI_COMM_WORLD);
        }


        fclose(fp);
        //}
        meshQualityEvaluation(submesh, id, OUTPUT_PATH);


    }

    if(id == 0) cout << "successful!!!" << endl;
    MPI_Finalize();

    return 0;
}

