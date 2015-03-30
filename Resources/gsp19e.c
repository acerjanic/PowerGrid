/*
 *  Gridding Image Reconstruction for Spiral k-space Acquisitions
 *                           (c) 1999-2003
 *  Copyright by Douglas C. Noll and the Regents of the University of Michigan 
 *                           (c) 1992-1995
 *  Copyright by Douglas C. Noll and the University of Pittsburgh 
 *
 *  Algorithm is based on O'Sullivan, IEEE Trans. Med. Imag., 4(4):200-207, 
 *    1985 as refined by Jackson, et al., Trans. Med. Imag., 10(2):473-478, 
 *    1991.  Sample density correction as described by Meyer, et al., 
 *    Magn. Reson. in Med., 28:202-213, 1992. Time-segmented inhomogeneity 
 *    correction as described by Noll, et al., IEEE Trans. Med. Imag., 
 *    10(4):629-637, 1991. Spiral gradient design described by Glover in 
 *    Magn. Res. Med. 42(2):412-5, 1999
 *
 *  John Pauly of Stanford University wrote code that served as the
 *    basis for this program.  Many changes occurred since then.
 *  Several utility routines came from Numerical Recipes in C,
 *    namely bessi0 and four1.
 *  rdbm.h is a GE header file necessary to extract info from raw data files
 *  getrttrajghg.c generates spiral gradients (derived from G. Glover's code)
 *  getrttrajvd.c generates spiral gradients with variable density 
 *    (derived for from C. Meyer's code)
 *  header.h specifies the ANALYZE file format - from AIR (R. Woods)
 *  bio.h,libbio.c does big/little endian byte swapping - from G. Hood
 *    at Pittsburgh Supercomputer Center 
 *  
 *  Versions:
 *     gpr2   - converted to spin-echo proj.
 *     gpr3   - window = 2, added apodization correction
 *     gpr6   - added other convolution windows as well
 *              option for sample density correction
 *     gpr8   - added k-b window
 *     gsp1   - 1st run at spiral recon
 *     gsp2   - fixed set up for multi-slices
 *     gsp3   - added multi-phase, reordered input parm list
 *     gsp4   - added over flag (=2 for 2x over, =1 for no over)
 *            - works, but not as well as Craig's version
 *	      - fixed delay term on readout
 *     gsp4hc - 1st run at homogeneity correction (with field map)
 *     gsp4hc2 - uses a previously generated field map
 *     gsp4hc3 - added fixviews subroutine for phase corrections
 *     gsp4hc4 - added section to read in list of files and uncompress
 *     gsp7  - put shifting, rotating, etc. into gridding subroutine
 *           - brought map making into same program (-m) option
 *           - defaut is no correction, (-h) option will do correction
 *           - compressed input is made an option (-c)
 *     gsp9  - 5.x version
 *     gsp10 - added parameters to rhuser variables in header
 *     gsp11 - made subsequent times through faster than first 
 *     gsp12 - added phase multiplier
 *           - turned off most output unless using the (-v) option
 *           - added image registration (-r)
 *     gsp14 - added multi-coil recon option
 *           - added (-2) option to double translations with (-r)
 *           - added (-o) option to specify starting image number
 *           - added (-m2,-h2) options to do linear term inhomogeneity corr.
 *           - added (-m3,-h3) options for both linear and time-seg corrs.
 *           - added (-l) option for shift correction in obliques
 *     gsp15 - realtime grad generation hooks from craig
 *           - fixed bug in (-l and -r) options
 *     gsp16 - LX (8.x) version
 *           - modified parameter list for glover pulse sequence
 *           - added hooks to getrttrajghg - Glover's grad code
 *     gsp16a - handling split data acq, fixed (-l) option for 8.x
 *     gsp16b - fixed field mapping, inhomogeneity corrections
 *     gsp17  - added output code for analyze format (-A)
 *              (uses Roger Woods header.h from AIR package)
 *            - added -f flag to flip (lr) to convert from Rad to Neuro coords
 *            - added -R flag to rotate images (SPM users requested this)
 *            - added printouts of more header info
 *     gsp18  - did some modest house cleaning...
 *            - autodetect, read LX header and data byte swapped
 *            - removed all stdin inputs, everything is now command line
 *            - fixed phase array
 *     gsp18a - fixed negative delays in -d flag
 *            - output slice location info with -l
 *     gsp18b - added -fx and -fy flags to do flips for analyze
 *     gsp18d - added -I flag for inteleaved ordering of slices
 *            - added -C flag to adjust output scaling
 *     gsp18e - added variable desnity gradients
 *            - added sensefact var
 *     gsp18f - added -U flag for UNFOLD recon
 *     gsp19a - changed disdaq handling - field map now a disdaq
 *            - added rev spiral 
 *     gsp19b - added -Q flag for out recon in in-out spiral
 *     gsp19c - modified gridding to remove post-grid dens compensation
 *     gsp19d - fixed (maybe) pixel shift problem for arbitrary rots/transposes
 *     gsp19e - added code to do L2 combination of in-out spirals
 *              
 */
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rdbm.h> /* LX header info */
#include <header.h> /* ANALYZE header info */
#include <bio.h> /* PSC code to do byte swapping */
/* LX image header stuff */
#define IM_CTR_R 136            /*CENTER R COORD OF PLANE IMAGE*/
#define IM_CTR_A 140            /*CENTER A COORD OF PLANE IMAGE*/
#define IM_CTR_S 144            /*CENTER S COORD OF PLANE IMAGE*/
#define IM_NORM_R 148           /*NORMAL R COORD*/
#define IM_NORM_A 152           /*NORMAL A COORD*/
#define IM_NORM_S 156           /*NORMAL S COORD*/
#define IM_TLHC_R 160           /*R COORD OF TOP LEFT HAND CORNER*/
#define IM_TLHC_A 164           /*A COORD OF TOP LEFT HAND CORNER*/
#define IM_TLHC_S 168           /*S COORD OF TOP LEFT HAND CORNER*/
#define IM_TRHC_R 172           /*R COORD OF TOP RIGHT HAND CORNER*/
#define IM_TRHC_A 176           /*A COORD OF TOP RIGHT HAND CORNER*/
#define IM_TRHC_S 180           /*S COORD OF TOP RIGHT HAND CORNER*/
#define IM_BRHC_R 184           /*R COORD OF BOTTOM RIGHT HAND CORNER*/
#define IM_BRHC_A 188           /*A COORD OF BOTTOM RIGHT HAND CORNER*/
#define IM_BRHC_S 192           /*S COORD OF BOTTOM RIGHT HAND CORNER*/
/* my defines */
#define MXNPR 16
#define MXNIM 256
#define MXNDAT 16384
#define MXNSL 64
#define NWEIGHTS 1024
#define SEGSIZE 3000 /* us */ 
#define OVER 2
#define FSIZE 3
#define EPS 0.0000001
#define GRIDL 1.5
#define GRIDB 6.7  /* Jackson's betas 1.5 = 6.7, 2.0 = 9.1, 2.5 = 11.5 */
#define USC (unsigned char *) 
#define PI    3.14159265358979323846
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* Raw Header stuff */
POOL_HEADER header;
RDB_HEADER_REC *h_rec = &header.rdb_hdr_rec;
EXAMDATATYPE *e_rec = (EXAMDATATYPE *)header.rdb_hdr_exam;
SERIESDATATYPE *s_rec = (SERIESDATATYPE *)header.rdb_hdr_series;
MRIMAGEDATATYPE *i_rec = (MRIMAGEDATATYPE *)header.rdb_hdr_image;
RDB_SLICE_INFO_ENTRY *a_rec = (RDB_SLICE_INFO_ENTRY *)header.rdb_hdr_data_acq_tab;
/* arrays and variables  */
float refbuf[MXNSL][MXNIM][MXNIM];
float reflbuf[MXNSL][3];
short int outbuf[MXNSL][MXNIM][MXNIM];
short int outbuf2[MXNSL][MXNIM][MXNIM];
short int ib1[MXNDAT][2];
float w1[2][MXNIM*OVER];
float prd[2][MXNPR][MXNDAT];
float im[2][MXNIM*OVER][MXNIM*OVER];
float sampim[MXNIM*OVER][MXNIM*OVER];
float grim[2][MXNIM*OVER][MXNIM*OVER];
float weight[NWEIGHTS+1];
float refim[MXNIM][MXNIM];
float refimmag[MXNIM][MXNIM];
float finalim[2][MXNIM][MXNIM];
float *t2k[2][MXNPR];
float ws[MXNIM*OVER][MXNIM*OVER];
float kdens[MXNDAT];
float vdfact[MXNDAT];
float regxs[40][2000],regys[40][2000],regrot[40][2000];
float refl[3];
char kroot[256];
char refroot[256];
char uncfile[256];
char regfile[256];
char com[256];
char datestr[12];
FILE *fk, *fi, *fo1, *fr;
int numrot=0, nim=64, samp1=0, samp1a, isl=0, iph=0, fsize=FSIZE;
int chop, samp_cor, file, nex;
int npr, ndat, ndatfr, nsamp, nslices, nsegs,exam,durmin,dursec;
int slnum,phnum,nphases,nph1,nphmult,rawnum,rotation,transpose;
int coilnum,ncoils,mcskip,imoffset=0;
int map_flg=0,hc_flg=0,com_flg=0,first_time=1,phs_flg=0,out_flg=0,loc_flg=0;
int allcoils_flg=0,reg_flg=0,lincor_flg=0,linmap_flg=0,svraw_flg=0;
int rev_flg = 0,anal_flg=0, flipx_flg=0, flipy_flg=0, nim_flg=0, unf_flg=0, unfact=1;
int sliceorder=1, rsp_flg=0, rsa_flg=0, L2_flg=0, norec_flg=0;
float lrshift=0.0,tbshift=0.0,zoomer=1.0;
float pt=0.0,magfact=0.0,phasefact=0.0;
float wind,gridl,gridb,maxk,slthick;
float factxx,factxy,factyx,factyy;
float factxxz,factxyz,factyxz,factyyz;
float pix_shifth,pix_shiftv;
float outscale = 512.0;
float sensefact=1.0;
double ts,gts,fsgcm,opfov;
int risetime, densamp, mapdda, dda, phnumout;
int gtype, opxres, ngap, concat, byteswap, densamp;
int stph, enph;
double slewrate, fast_rec_off, mapdel, samptime;
int fast_rec_lpf;
/* functions */
double rint();
float kaiser();
void fft();
void loc_calc();

main(argc,argv)
int argc;
char **argv;
{
  int i,j,k,stsl,ensl,segnum,midsamp,sampwind;
  int tmpslice,pp,ss;
  float jj;
  char *fnb;

  if (argc< 2) {
fprintf(stderr,"Usage: %s [OPTIONS] rawfile1 rawfile2...\n",argv[0]);
fprintf(stderr,"\nRecon Options\n");
fprintf(stderr,"-h   Do inhomogeneity correction (segmented)\n");
fprintf(stderr,"-h2  Do inhomogeneity correction (linear only)\n");
fprintf(stderr,"-h3  Do inhomogeneity correction (-h and -h2)\n");
fprintf(stderr,"-m   Generate field maps (for use with -h)\n");
fprintf(stderr,"-m2  Generate field maps (for use with -h2)\n");
fprintf(stderr,"-m3  Generate field maps (for use with -h3)\n");
fprintf(stderr,"-F # Filter size for smoothing fieldmap (default %d)\n",FSIZE);
fprintf(stderr,"-c   Compressed raw files\n");
fprintf(stderr,"-M # Navigator correction for multishot (magnitude)\n");
fprintf(stderr,"-P # Navigator correction for multishot (phase)\n");
fprintf(stderr,"-O # Off-resonance demodulation (cycles over readout)\n");
fprintf(stderr,"-n # Recon image size (default set by pulse sequence)\n");
fprintf(stderr,"-d # Sample delay term (for calibration)\n");
fprintf(stderr,"-r F Does image registration according to file name \"F\"\n");
fprintf(stderr,"-U   UNFOLD recon\n");
fprintf(stderr,"-Q   Out spiral data in in-out spiral acquisition\n");
fprintf(stderr,"\nOutput Options\n");
fprintf(stderr,"-V   Read header only (no recon)\n");
fprintf(stderr,"-v   Verbose\n");
fprintf(stderr,"-l   Shift images to actual Rx'ed location\n");
fprintf(stderr,"-A   Output Analyze format\n");
fprintf(stderr,"-R1  Rotation image 90 deg.\n");
fprintf(stderr,"-R2  Rotation image 180 deg. (also -R)\n");
fprintf(stderr,"-R3  Rotation image 270 deg.\n");
fprintf(stderr,"-fx  Flip images in x-direction\n");
fprintf(stderr,"-fy  Flip images in y-direction\n");
fprintf(stderr,"-q   Reverse slice ordering\n");
fprintf(stderr,"-x # Shift images in x-direction by # pixels\n");
fprintf(stderr,"-y # Shift images in y-direction by # pixels\n");
fprintf(stderr,"-z # Zoom factor (default 1.0)\n");
fprintf(stderr,"-p   Output phase images (radians*1000)\n");
fprintf(stderr,"-a   Output images for all coils in multicoil recon\n");
fprintf(stderr,"-L   Do L2 norm with existing images (for spiral in/out)\n");
fprintf(stderr,"-t # Recon temporal frame # only\n");
fprintf(stderr,"-s # Recon slice # only (incompatible with -A)\n");
fprintf(stderr,"-o # Offset output image numbers (adds # to temp frame)\n");
fprintf(stderr,"-S   Output mag k-space data (for diags)\n");
fprintf(stderr,"-I   For interleaved slice ordering\n");
fprintf(stderr,"-C # To adjust output scaling (default = %f)\n",outscale);
    return(0);
  }

  file = 1;
  while (argv[file][0] == '-')
    switch (argv[file++][1]) { 
      case 'h':                          /* inhomogeneity correction */
        switch (argv[file-1][2]) {
	  case '2': lincor_flg = 1; break;  /* const, linear terms */
	  case '3': lincor_flg = 1; hc_flg = 1; break; /* both ways */
          default: hc_flg = 1; break;    /* time-segmented only */
	} break;
      case 'm':                          /* correction maps for hom. corr. */
        switch (argv[file-1][2]) {
	  case '2': linmap_flg = 1; break;  /* const, linear terms */
	  case '3': linmap_flg = 1; map_flg = 1; break; /* both ways */
          default: map_flg = 1; break;   /* gen. field map */
	} break;
      case 'R':                          /* extra rotations */
        switch (argv[file-1][2]) {
	  case '1': numrot = 1; break;  /* 90 deg.*/
	  case '2': numrot = 2; break;   /* 180 deg. */
	  case '3': numrot = 3; break;   /* 270 deg. */
          default: numrot = 2; break;    /* defaults to 180 */
	} break;
      case 'c': com_flg = 1; break;      /* compressed input */
      case 'p': phs_flg = 1; break;      /* output phase images */
      case 'v': out_flg = 1; break;      /* verbose - output status info */
      case 'V': out_flg = 1; norec_flg = 1;break;  /* verbose - output status info */
      case 'a': allcoils_flg = 1; break; /* output images for all coils */
      case 'L': L2_flg = 1; break;       /* output L2 combintation of images */
      case 'r': reg_flg = 1; strcpy(regfile,argv[file++]); break;
					 /* image registration */
      case 'o': imoffset = atoi(argv[file++]); break;  
					 /* offset in output image numbers */
      case 'n': nim_flg = 1; nim = atoi(argv[file++]); break;  
      case 't': iph = atoi(argv[file++]); break;  
      case 's': isl = atoi(argv[file++]); break;  
      case 'd': samp1 = atoi(argv[file++]); break;  
      case 'F': fsize = atoi(argv[file++]); break;  
      case 'x': lrshift = atof(argv[file++]); break;  
      case 'y': tbshift = atof(argv[file++]); break;  
      case 'z': zoomer = atof(argv[file++]); break;  
      case 'C': outscale = atof(argv[file++]); break;  
      case 'O': pt = atof(argv[file++]); break;  
      case 'P': phasefact = atof(argv[file++]); break;  
      case 'S': svraw_flg = 1; break;    /* save mag raw data into images */
      case 'l': loc_flg = 1; break;      /* ture locations, shift images */
      case 'q': rev_flg = 1; break;      /* reverse slice ordering */
      case 'A': anal_flg = 1; break;     /* ANALYZE output format */
      case 'U': unf_flg = 1; break;      /* UNFOLD recon */
      case 'Q': rsa_flg = 1; break;      /* Out flag in in-out spiral */
      case 'I': sliceorder = 0; break;   /* 0=interleaved, 1=sequential */
      case 'f': 
        switch (argv[file-1][2]) {
	  case 'x': flipx_flg = 1; break;     /* flip x axis at output */
	  case 'y': flipy_flg = 1; break;     /* flip y axis at output */
          default:  flipx_flg = 1; break;     /* flip x axis (default) */
	} break;
    
      default: fprintf(stderr,"%s: Did not recognize %s",argv[0],argv[file-1]);
    }

/* printf("phasefact = %f\n",phasefact); */

 for (rawnum = 0; rawnum < (argc - file); rawnum++)
 {
  fnb = argv[rawnum+file];
  if (com_flg) {
    sprintf(uncfile,"%s.uncp", fnb);
    if (out_flg) printf("Processing raw file %s\n", fnb);
    sprintf(com,"uncompress -c %s > %s", fnb, uncfile);
    system(com);
  }
  else
    sprintf(uncfile,"%s", fnb);

  /* get some parameters from header of first file */
  if (rawnum == 0) {
    get_hdr_info();
    cal_map();
    if (map_flg) phasefact = 0.0;
    if (reg_flg) {
      if (!(fr=fopen(regfile,"r"))) {
          fprintf(stderr,"Can't open %s!\n",regfile); exit(0);
      }
      while (fscanf(fr,"%d %d",&pp,&ss)==2) {
        fscanf(fr,"%f %f %f %f",&regxs[ss][pp],&regys[ss][pp],&regrot[ss][pp],&jj);
	/* printf("%d %d %f %f %f %f\n",pp,ss,regxs[ss][pp],regys[ss][pp],regrot[ss][pp],jj); */
      }
      fclose(fr); 
    }
  }

  if (norec_flg) break;

  for (coilnum = 0; coilnum < ncoils; coilnum ++){

  if (map_flg || linmap_flg) {
    /* do field map */
    stph = 0; enph = 2; iph = 0; unf_flg = 0;
    for (slnum = 0; slnum < nslices; slnum++) {
      for(phnum = stph; phnum < enph; phnum++) {
      if (out_flg) { 
	printf("Co %d, Ph %2d, Sl %2d: .",coilnum+1,phnum+1,slnum+1);
        fflush(stdout); }
      if(!load_projs()) break;
      refocus();
      if (out_flg) { printf("."); fflush(stdout); }
      load_ft1();
      if (out_flg) { printf("."); fflush(stdout); }
      fermi_filt1(); 
      fermi_filt2(); 
      if (out_flg) { printf("."); fflush(stdout); }
      if (hc_flg) {
	/* segmented hc recon */
        nsegs = floor(1.0*ndat*samptime/SEGSIZE);
        for (segnum = 0; segnum <= nsegs; segnum++) {
          if (out_flg) { printf("%d.",segnum+1); fflush(stdout); }
	  sampwind = SEGSIZE/samptime;
	  midsamp = segnum*sampwind;
	  wind_dat(midsamp,sampwind);
          ft_image();
	  make_final(segnum);
          /* write_image(fnb,slnum+1,segnum+1); */
        }
      }
      else {
	/* plain old recon */
        copy_dat();
        ft_image();
      }
      if (out_flg) { printf("."); fflush(stdout); }
	if (phnum == 0) {
	  make_phase1_filt();
	}
	else {
	  make_phase2_filt(); 
	  if (ncoils > 1) {
	    write_ref_mc(slnum,coilnum);
	    if (coilnum == (ncoils-1)) {
	      write_ref_mcf(slnum);
	      write_ref(slnum);
	    }
	  } else 
	    write_ref(slnum);
        }
      if (out_flg) { printf("\n"); }
      }
    }
 } else {
  /* do recon */
  if (unf_flg) unfact = npr; else unfact = 1;
  if (iph == 0) {  stph = (mapdda*unfact); enph = nphases*unfact; }
  else {  stph = iph-1+mapdda; enph = iph+mapdda; }
  if (isl == 0) {  stsl = 0; ensl = nslices; }
  else {  stsl = isl-1; ensl = isl; }
  if (hc_flg || lincor_flg) load_ref_all(); 
   for(phnum = stph; phnum < enph; phnum++) {
    phnumout = phnum+1-mapdda*unfact;
    for (slnum = stsl; slnum < ensl; slnum++) {
      if (hc_flg || lincor_flg) load_ref(slnum); 
      if (out_flg) { 
	printf("Co %d, Ph %2d, Sl %2d: .",coilnum+1,phnumout,slnum+1);
        fflush(stdout); }
      if(!load_projs()) break;
      refocus();
      if (npr > 1) fixviews(); 
      if (out_flg) { printf("."); fflush(stdout); }
      if (first_time || reg_flg || lincor_flg) {
        load_ft1();
        if (out_flg) { printf("."); fflush(stdout); }
        fermi_filt1(); 
        fermi_filt2(); 
        first_time = 0;
      }
      else {
        load_ft2();
        if (out_flg) { printf("."); fflush(stdout); }
        fermi_filt2(); 
      }
      if (svraw_flg) write_raw(); 
      if (out_flg) { printf("."); fflush(stdout); }
      if (hc_flg) {
	/* segmented hc recon */
        nsegs = floor(1.0*ndat*samptime/SEGSIZE);
        for (segnum = 0; segnum <= nsegs; segnum++) {
          if (out_flg) { printf("%d.",segnum+1); fflush(stdout); }
	  sampwind = SEGSIZE/samptime;
	  midsamp = segnum*sampwind;
	  wind_dat(midsamp,sampwind);
          ft_image();
	  make_final(segnum);
          /* write_image(fnb,slnum+1,segnum+1); */
        }
      }
      else {
	/* plain old recon */
        copy_dat();
        ft_image();
      }
      if (flipx_flg) flipx_image();
      if (flipy_flg) flipy_image();
      if (out_flg) { printf("."); fflush(stdout); }
      if (anal_flg)
        write_anal(slnum+1);
      else {
        /* write_image(slnum+1,phnumout,rawnum); */
        if (allcoils_flg) write_image(coilnum+1,slnum+1,phnumout,rawnum);
        write_imagemc(coilnum+1,slnum+1,phnumout,rawnum);
        if (phs_flg) write_image_phs(slnum+1,phnumout,rawnum);
      } 
      if (out_flg) { printf("\n"); }
    } /* slices loop */
    if (anal_flg)
      write_anal_all(coilnum+1,phnumout,rawnum);
   } /* phases loop */
   } /* if domap */
  } /* coils loop */
  if (com_flg) {
    sprintf(com,"rm -f %s", uncfile);
    system(com);
  }
 } /* raw file loop */
  if (out_flg) { printf("Done!!\n"); }
} /* main */

/* 8x header info from first file on list */
get_hdr_info()
{
    int i,k,n;
    int tmpi;
    float tmp1, tmp2;
    char *imptr;
    short int *shptr;
    long *intptr;
    float *flptr, *flptr2;
    float headflt();
    float s_i_offset,tempr2,tempr3;
    float p1[3],p2[3],p3[3],c[3],n2[3],n3[3];

    if (!(fi=fopen(uncfile,"r"))) {
      fprintf(stderr,"Can't open %s!\n",uncfile);
      return(0);
    }

    /* header */
    InitBIO();
    bio_big_endian_input = 1;
    if (!bio_big_endian_machine)
      if (out_flg) printf("Byte order swapping required (bio)\n\n");
    fread(&header,sizeof(header),1,fi);
    /* if ( (*h_rec).rdb_hdr_npasses > 255 ) 
    {
      byteswap = 1;
      if (out_flg) printf("Byte order swapping required\n\n");
    } else byteswap = 0; */
    nph1 = 1;
    nphmult = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user1)); 
    if (nphmult < 1) nphmult = 1;
    nphases = nph1*nphmult;
    if ((map_flg || linmap_flg)  && (nphases < 2)) {
      fprintf(stderr,"Illegal mapping file %s!",uncfile);
      exit(0);
    }
    if (map_flg || linmap_flg) nphases = 2;
    nslices = BRdInt16(USC&(*h_rec).rdb_hdr_nslices)/nph1;
    ndatfr = BRdInt16(USC&(*h_rec).rdb_hdr_frame_size);
    npr = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user4)); 
    dda = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user16)); 
    mapdda = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user19)); 
    chop = 1;
    nex = 1;
    densamp = 5;
    /* glover cv's */
    gtype = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user5)); /* future use */
    opxres = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user3)); 
    slewrate = BRdFloat32(USC&(*h_rec).rdb_hdr_user7);
    /* var dens cv's */
    densamp = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user17)); 
    sensefact = BRdFloat32(USC&(*h_rec).rdb_hdr_user18); 
    if (sensefact < 0.01) sensefact = 1.0;
    /* concat read */
    ngap = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user10)); 
    concat = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user11)); 
    tmpi = (int) rint(BRdFloat32(USC&(*h_rec).rdb_hdr_user21)); 
    /* reverse spiral flags
     * tmpi   rsp_flg   rsa_flg  concat
     *  0       0          0       0/1   standard (forward spiral) pulse seq
     *  1       1          0        0    reverse only 
     *  2       1          0        1    for/reverse - reconstruct rev 
     *  2       0          1        1    for/reverse - reconstruct for
     */
    if (tmpi > 0) {     
      ndat = ndatfr;    
      if (tmpi == 1)  { 
        rsp_flg = 1;
        rsa_flg = 0;
      } else            
	  rsp_flg = 1-rsa_flg;
    } else {            
      rsp_flg=0;
      rsa_flg=0;
      ndat = ndatfr*(concat+1);
    }
    /* fast rec cutoff kHz */
    fast_rec_off = BRdFloat32(USC&(*h_rec).rdb_hdr_user12); 
    fast_rec_lpf = BRdInt32(USC&(*h_rec).rdb_hdr_fast_rec);
    mapdel = BRdFloat32(USC&(*h_rec).rdb_hdr_user15); /* field map offset (us) */
    samptime = BRdFloat32(USC&(*h_rec).rdb_hdr_user13); /* in (us) */
    fsgcm = BRdFloat32(USC&(*h_rec).rdb_hdr_user6);
    risetime = (int) rint(fsgcm*10000.0/slewrate);
    opfov = BRdFloat32(USC&(*h_rec).rdb_hdr_user0);
    ts = BRdFloat32(USC&(*h_rec).rdb_hdr_user13)*1e-6;  /* in sec */
    gts = 4*1e-6;  /* 4 us gradient spacing in sec */
    if (!nim_flg) nim = BRdInt16(USC&(*h_rec).rdb_hdr_im_size);
    ncoils = BRdInt16(USC&(*h_rec).rdb_hdr_dab[0].stop_rcv)-BRdInt16(USC&(*h_rec).rdb_hdr_dab[0].start_rcv)+1;
    mcskip = BRdInt32(USC&(*h_rec).rdb_hdr_raw_pass_size)/ncoils;
    imptr = (char *)i_rec;
    shptr = (short int *)i_rec;
    intptr = (long *)i_rec;
    flptr = (float *)i_rec;
    flptr2 = (float *)a_rec;
    exam = BRdInt16(USC(shptr+4));
    slthick = BRdFloat32(USC(flptr+7));
    strcpy(datestr,(*h_rec).rdb_hdr_scan_date);
    datestr[2] = '_';
    datestr[5] = '_';
    durmin = floor(BRdFloat32(USC(flptr+6))/60000000.0);
    dursec = ceil(BRdFloat32(USC(flptr+6))/1000000.0 - 60.0*durmin);
    rotation = BRdInt16(USC&(*h_rec).rdb_hdr_rotation); 
    transpose = BRdInt16(USC&(*h_rec).rdb_hdr_transpose); 

    /* for (i=0; i<256; i++)
      printf(" %d %d %d %d %f\n",i,*(shptr+2*i),*(shptr+2*i+1),*(intptr+i),*(flptr+i));  */
    if (out_flg) {
      printf("recon image matrix = %d\n",nim);
      printf("number of time points = %d (%d*%d)\n",nphases,nph1,nphmult);
      printf("number of slices = %d\n",nslices);
      printf("number of samples in each acq frame = %d\n",ndatfr);
      printf("number of samples in readout = %d\n",ndat);
      printf("number of disdaqs = %d + %d (map)\n",dda,mapdda);
      printf("rev spiral flag = %d\n",tmpi);
      printf("field of view (cm) = %.1f\n",opfov);
      printf("slice thickness (mm) = %.1f\n",slthick);
      printf("TR (ms) = %.1f\n",0.001*(BRdInt32(USC(intptr+50))));
      printf("TE (ms) = %.1f\n",0.001*(BRdInt32(USC(intptr+52))));
      printf("flip angle = %d\n",BRdInt16(USC(shptr+130)));
      printf("number of spirals = %d\n",npr);
      printf("time of scan = %s\n",(*h_rec).rdb_hdr_scan_time);  
      printf("date of scan = %s\n",datestr);
      printf("total scan duratios (mm:ss) = %2d:%02d\n",durmin,dursec);
      printf("pulse sequence name = %s\n",(imptr+320));
      printf("exam number = %d\n",exam);
      printf("size of header = %d\n", sizeof(header));
      printf("nphmult = %d\n",nphmult); 
      printf("nph1 = %d\n",nph1); 
      printf("ncoils = %d\n",ncoils); 
      printf("mcskip = %d\n",mcskip); 
      printf("fsgcm = %.1f, slewrate = %.1f, risetime = %d, opfov = %.1f\n", fsgcm, slewrate, risetime, opfov); 
      printf("sensefact = %f\n",sensefact);
      printf("concat = %d, ngap = %d\n",concat, ngap); 
      printf("opxres = %d\n",opxres);
      printf("mapdel = %f\n",mapdel);
      printf("fast_rec_lpf = %d, fast_rec_off = %f\n",fast_rec_lpf,fast_rec_off);
      printf("ts = %f (us), gts = %f (us)\n",ts*1e6,gts*1e6);
      printf("rot = %d\n",rotation);
      printf("trans = %d\n",transpose);
      printf("prescan r1 = %d\n",BRdInt32(USC&(*h_rec).rdb_hdr_ps_aps_r1)); 
      printf("prescan r2 = %d\n",BRdInt32(USC&(*h_rec).rdb_hdr_ps_aps_r2)); 
      printf("prescan tg = %d\n",BRdInt32(USC&(*h_rec).rdb_hdr_ps_aps_tg)); 
      printf("prescan freq = %d\n",BRdInt32(USC&(*h_rec).rdb_hdr_ps_aps_freq));
      printf("raw file size = %d\n",BRdInt32(USC&(*h_rec).rdb_hdr_raw_pass_size)+sizeof(header));
    }
    /* all of the old image location stuff (-l) appears to be broken in 8.x */
    /* and has been replaced by info passed from pulse sequence  */
    if (loc_flg) {
      pix_shifth = -BRdFloat32(USC&(*h_rec).rdb_hdr_yoff)/BRdInt16(USC&(*h_rec).rdb_hdr_im_size);
      pix_shiftv = BRdFloat32(USC&(*h_rec).rdb_hdr_xoff)/BRdInt16(USC&(*h_rec).rdb_hdr_im_size);
      if (out_flg) 
        printf("shifth = %f, shiftv = %f (in mm)\n",pix_shifth*opfov*10,pix_shiftv*opfov*10);
    }
    else
      pix_shifth = pix_shifth = 0.0;

      /*
      printf("corner locations for image 1:\n");
      printf(" TL:  %f, %f, %f\n", headflt(imptr + IM_TLHC_R), headflt(imptr + IM_TLHC_A), headflt(imptr + IM_TLHC_S));
      printf(" TR:  %f, %f, %f\n", headflt(imptr + IM_TRHC_R), headflt(imptr + IM_TRHC_A), headflt(imptr + IM_TRHC_S));
      printf(" BR:  %f, %f, %f\n", headflt(imptr + IM_BRHC_R), headflt(imptr + IM_BRHC_A), headflt(imptr + IM_BRHC_S)); */
    
    for (i = 0; i < 3; i++)
    {
        p1[i] = BRdFloat32(USC( &(*(a_rec)).gw_point1[i]));
        p2[i] = BRdFloat32(USC( &(*(a_rec)).gw_point2[i]));
        p3[i] = BRdFloat32(USC( &(*(a_rec)).gw_point3[i]));
    }
    loc_calc(p1,p2,p3,rotation,transpose);
    s_i_offset = BRdFloat32(USC(flptr+42)) - p1[2];
    /* calculate in-plane normal vectors n2,n3 */
    /* old method
    c[0] = BRdFloat32(USC(flptr+34));
    c[1] = BRdFloat32(USC(flptr+35));
    c[2] = BRdFloat32(USC(flptr+36));
    c[2] -= s_i_offset;
    tempr2 = tempr3 = 0.;
    for (i = 0; i < 3; i++)
    {
      n2[i] = p2[i] - p1[i];
      tempr2 += n2[i]*n2[i];
      n3[i] = p2[i] - p3[i];
      tempr3 += n3[i]*n3[i];
    }
    tempr2 = sqrt(tempr2);
    tempr3 = sqrt(tempr3);
    for (i = 0; i < 3; i++)
    {
      n2[i] /= tempr2; 
      n3[i] /= tempr3; 
    } 
    printf(" TL:  %f, %f, %f\n",p1[0],p1[1],p1[2]+s_i_offset); 
    printf(" TR:  %f, %f, %f\n",p2[0],p2[1],p2[2]+s_i_offset); 
    printf(" BR:  %f, %f, %f\n",p3[0],p3[1],p3[2]+s_i_offset); 
    printf(" C:  %f, %f, %f\n",c[0],c[1],c[2]+s_i_offset); 
    printf(" N2:  %f, %f, %f\n",n2[0],n2[1],n2[2]); 
    printf(" N3:  %f, %f, %f\n",n3[0],n3[1],n3[2]); 
    printf("tempr2 %f, tempr3 %f \n",tempr2,tempr3); */

    if (loc_flg) {
      /* old method 
      pix_shifth = -(n2[0]*c[0] +  n2[1]*c[1] +  n2[2]*c[2])/tempr2;
      pix_shiftv = (n3[0]*c[0] +  n3[1]*c[1] +  n3[2]*c[2])/tempr2;
      printf("  shifth = %f, shiftv = %f (in mm)\n",pix_shifth*tempr2,pix_shiftv*tempr2);
      */

      if (out_flg)
      for (n = 0; n < nslices; n++)
      {
        for (i = 0; i < 3; i++)
        {
          p1[i] = BRdFloat32(USC( &(*(a_rec+n)).gw_point1[i]));
          p2[i] = BRdFloat32(USC( &(*(a_rec+n)).gw_point2[i]));
          p3[i] = BRdFloat32(USC( &(*(a_rec+n)).gw_point3[i]));
        }

        loc_calc(p1,p2,p3,rotation,transpose);
        printf("slice %d    R/L       A/P       S/I\n",n+1);
        printf(" TL:  %f, %f, %f\n",p1[0],p1[1],p1[2]+s_i_offset);
        printf(" TR:  %f, %f, %f\n",p2[0],p2[1],p2[2]+s_i_offset);
        printf(" BR:  %f, %f, %f\n",p3[0],p3[1],p3[2]+s_i_offset);
      }
   }
   else
      pix_shifth = pix_shifth = 0.0;


    /* old 5.x stuff deleted */

    /* perform transposing and rotating */
    factxx = -1.0; factxy = 0.0; factyx = 0.0; factyy = -1.0; 
    if (transpose) {
      tmp1 = factyx; tmp2 = factyy;
      factyx = factxx; factyy = factxy;
      factxx = tmp1; factxy = tmp2;
    }
    for (k = 0; k < rotation; k++) {
      tmp1 = factyx; tmp2 = factyy;
      factyx = -factxx; factyy = -factxy;
      factxx = tmp1; factxy = tmp2;
    }

    /* do rotations on pixel shifts (this still may not be correct for all cases) */
    tmp1 = factxx*pix_shifth + factxy*pix_shiftv;
    tmp2 = factyx*pix_shifth + factyy*pix_shiftv;
    /* printf("shifth = %f, shiftv = %f; new shifth = %f, shiftv = %f\n",pix_shifth,pix_shiftv,tmp1,tmp2); */
    pix_shifth = tmp1;
    pix_shiftv = tmp2;

    fclose(fi); 

    /* just temp
    if (!(fi=fopen("header","w"))) {
      fprintf(stderr,"Can't open %s!\n",uncfile);
      return(0);
    }

    (*h_rec).rdb_hdr_nslices = 12;

    fwrite(&header,sizeof(header),1,fi);
    fclose(fi);  */ 
}

cal_map()
{
  int i,n,j,k;
  int res;
  int getrttrajghg();
  int getrttrajvd();
  float tmp1,tmp2;
  double kx,ky,kxo,kyo,ggx,ggy,ggm,kksr,kksi;

  /* set-up zoom and rotation parameters */
  factxxz = 1.0/zoomer; factxyz = 0.0; factyxz = 0.0; factyyz = 1.0/zoomer; 
  for (k = 0; k < (numrot); k++) {
    tmp1 = factyxz; tmp2 = factyyz;
    factyxz = -factxxz; factyyz = -factxyz;
    factxxz = tmp1; factxyz = tmp2;
  }

  /* kaiser-bessel functions - initialize terms */
  gridl = GRIDL; gridb = GRIDB; wind = gridl;

  refl[0] = refl[1] = refl[2] = 0.0;

  if (gtype==0) {
    res = getrttrajghg(opxres, npr, ts, gts, fsgcm, opfov, slewrate, 
		gts*21000, gtype, t2k[0], t2k[1]);
    for (i=0; i < res; i++)
	    vdfact[i]=1;
  } else {
    res = getrttrajvd2(opxres, npr, densamp, sensefact, ts, gts, fsgcm, opfov, slewrate, gts*21000, gtype, t2k[0], t2k[1],vdfact);
  }
  
  if (out_flg) 
  {
    /* printf("ndat = %d, npr = %d, ts = %f, gts = %f, fsgcm = %f, opfov = %f, risetime = %d\n",ndat, npr, ts, gts, fsgcm, opfov, risetime); */
    printf("acq. matrix = %d, nom. resolution = %f mm\n",opxres,10*opfov/opxres); 
    printf("k-space points = %d\n",res); 
  }
  nsamp = ndat;

  /* sample density correction */
  /* Using Craig's formula from MRM, 28, 202-213 (1992) */
  kdens[0] = 0.0; 
  for (i=1; i < nsamp; i++)
  {
      /* might wish to add correction for delay */
      kx = t2k[0][0][i];
      ky = t2k[1][0][i];
      kxo = t2k[0][0][i-1];
      kyo = t2k[1][0][i-1];
      kksr = kx*kxo + ky*kyo;
      kksi = ky*kxo - kx*kyo;
      /* printf("%d: %f %f %f %f %f %f\n",i,kx,ky,kxo,kyo,kksr,kksi); */
      ggx = kx - kxo;
      ggy = ky - kyo;
      if (((ggm = hypot(ggx,ggy)) > 0.) && (hypot(kksr,kksi) > 0.) ) {
	tmp1 = ggm * fabs(sin(atan2(ggy,ggx) - atan2(ky,kx))); 
        /* tmp2 = fabs((hypot(kx,ky) - hypot(kxo,kyo))/atan2(kksi,kksr)); */
	kdens[i] = tmp1*vdfact[i]/310; /* 310 is a correction factor */
      } else
        kdens[i] = 0.;
     /* printf("%f %f %f %f %f\n", kx, ky,kdens[i],tmp2,vdfact[i]); */
    /* components to density correction:
      ggm = gradient strength; linear velocity
      sin(atan2(ggy,ggx) - atan2(ky,kx)) =  tangential outward of velocity 
             These first two came from the Meyer paper.
      vdfact - density correction term
      (hypot(kx,ky) - hypot(kxo,kyo))/atan2(kksi,kksr) = radial density
	     (distance outward)/(anglar distance along arc)
	     This term is different from the previous term in that this
	     measure true line-line distances in the radial direction.  
	     It was added for variable density trajectories (except for the
	     the first 2 or 3 points, it is constant for uniform density
	     spiral trajectories. 
    */
  }
  for (i=nsamp; i < ndat; i++) kdens[i] = 0.0; 

  if (samp1 < 0)
  {
    samp1a = -samp1;
    samp1 = 0;
  }
}

load_projs()
{
  int i, j, k, n, m, chopper=1;
  FILE *jj;
  int nodd;
  int coiloff, nproj, offs, offs2;
  int phnum2, pr2;

  nodd = (nphmult*npr*(concat+1))%2;

  /* openfile and head header */
  if (((slnum == 0) || (isl != 0)) && (phnum == stph) && (coilnum == 0))
  {
      if (!(fi=fopen(uncfile,"r"))) {
        fprintf(stderr,"Can't open %s!\n",uncfile);
        return(0);
      }

      /* header */
      fread(&header,sizeof(header),1,fi);
  }
 
  /* for UNFOLD */
  if ((unf_flg) && (unfact > 1)) {
    pr2 = phnum%npr; 
    phnum2 = (phnum-pr2)/npr;
  } else {
    pr2 = npr-1;
    phnum2 = phnum;
  }

  /* build pointer location - echo 1 */
  /* assumes nph1 = 1 */
  coiloff = mcskip*coilnum;
  nproj = slnum*((concat+1)*nphmult*npr+1+nodd) + (concat+1)*phnum2*npr + 1;
  offs = sizeof(header) + nproj*ndatfr*4 + coiloff;
  /* printf("slnum = %d, phnum = %d, coilnum = %d, nproj = %d, offs = %d\n",slnum,phnum,coilnum,nproj,offs);  */

  /* read 1st half of data */
  fseek(fi, offs, SEEK_SET);
  for (n=0; n<npr; n++) {
    FRdInt16Array(fi,ib1[0],2*ndatfr);
    if (bio_error) { printf("bio_error\n"); return(0); }
    if ((unf_flg) && (unfact > 1) && (n != pr2)) 
    {
        for (j=0; j<ndatfr; j++) {
          prd[0][n][j] = 0.;
          prd[1][n][j] = 0.;
        }
    }
    else
    {
      if ((rsp_flg==0) && (rsa_flg==0)) /* std. forward spiral */
        for (j=0; j<ndatfr; j++) {
          prd[0][n][j] = ib1[j][0]*chopper;
          prd[1][n][j] = ib1[j][1]*chopper;
        }
      else if (rsa_flg==0) /* reverse spiral */
        for (j=0; j<ndatfr; j++) {
          prd[0][n][ndatfr-j-1] = ib1[j][0]*chopper;
          prd[1][n][ndatfr-j-1] = ib1[j][1]*chopper;
        }
    }
    if (concat) {
      FRdInt16Array(fi,ib1[0],2*ndatfr);
      if (bio_error) { printf("bio_error\n"); return(0); }
      if ((unf_flg) && (unfact > 1) && (n != pr2)) 
        for (j=0; j<ndatfr; j++) {
          prd[0][n][j+ndatfr+ngap] = 0.;
          prd[1][n][j+ndatfr+ngap] = 0.;
        }
      else
      {
        if ((rsp_flg==0) && (rsa_flg==0)) /* std. forward spiral */
          for (j=0; j<ndatfr; j++) {
            prd[0][n][j+ndatfr+ngap] = ib1[j][0]*chopper;
            prd[1][n][j+ndatfr+ngap] = ib1[j][1]*chopper;
          }
        else if (rsa_flg==1) /* forward spiral from in-out */
          for (j=0; j<ndatfr; j++) {
            prd[0][n][j] = ib1[j][0]*chopper;
            prd[1][n][j] = ib1[j][1]*chopper;
          }
      }
    }

  /*
    if ((i=fread(ib1,4,ndatfr,fi))!=ndatfr) return(0);
    for (j=0; j<i; j++) {
      prd[0][n][j] = ib1[j][0]*chopper;
      prd[1][n][j] = ib1[j][1]*chopper;
    }
    if (concat) {
      if ((i=fread(ib1,4,ndatfr,fi))!=ndatfr) return(0);
      for (j=0; j<i; j++) {
        prd[0][n][j+ndatfr+ngap] = ib1[j][0]*chopper;
        prd[1][n][j+ndatfr+ngap] = ib1[j][1]*chopper;
      }
    }
    */
    chopper *= chop;
  }

  if (((phnum2 == (enph - 1) && (pr2 == (npr-1))) || (iph != 0)) && ((slnum == (nslices -1)) || (isl != 0)) && (coilnum == (ncoils - 1)))
    fclose(fi); 

  return(1);
}

refocus()
{
  int i, j, j2;
  float tr, ti,cs,sn,pt2,freqpt;

  if (rsp_flg==0)
    freqpt = ((fast_rec_off*1e3*ts) + pt/ndat +refl[0]/(2.0*PI));
  else
    freqpt = -((fast_rec_off*1e3*ts) + pt/ndat +refl[0]/(2.0*PI));
  for (j=0; j<ndatfr; j++) {
    cs = cos(-2.0*PI*freqpt*j);
    sn = sin(-2.0*PI*freqpt*j);
    for (i=0; i<npr; i++) {
      tr = prd[0][i][j]; ti = prd[1][i][j];
      prd[0][i][j] = tr*cs - ti*sn;
      prd[1][i][j] = tr*sn + ti*cs;
      /* printf("%f %f\n", prd[0][i][j],prd[1][i][j]);  */
    }
  }
  if ((concat) && ((rsp_flg+rsa_flg)==0)) {
    for (j=0; j<ndatfr; j++) {
    j2 = j+ndatfr+ngap;
    cs = cos(-2.0*PI*freqpt*j2);
    sn = sin(-2.0*PI*freqpt*j2);
    for (i=0; i<npr; i++) {
      tr = prd[0][i][j2]; ti = prd[1][i][j2];
      prd[0][i][j2] = tr*cs - ti*sn;
      prd[1][i][j2] = tr*sn + ti*cs;
      /* printf("%f %f\n", prd[0][i][j2],prd[1][i][j2]); */
    }
    }
    /* fill in missing points in views for concat readouts */
    for (i=0; i<npr; i++) {
    for (j=0; j<ngap; j++) {
      prd[0][i][ndatfr+j] = 0.5*(prd[0][i][ndatfr-1] + prd[0][i][ndatfr+ngap]);
      prd[1][i][ndatfr+j] = 0.5*(prd[1][i][ndatfr-1] + prd[1][i][ndatfr+ngap]); 
    }
    }
    /* for (j=ndatfr-5; j<ndatfr+5; j++) 
          printf("%f %f\n",prd[0][0][j],prd[1][0][j]); */
  }
  if (samp1a > 0)
  {
    for (j=ndat-samp1a-1; j>=0; j--) 
      for (i=0; i<npr; i++) {
        prd[0][i][j+samp1a] = prd[0][i][j];
        prd[1][i][j+samp1a] = prd[1][i][j];
    }
    for (j=1; j<samp1a; j++) 
      for (i=0; i<npr; i++) {
        prd[0][i][j] = prd[0][i][0];
        prd[1][i][j] = prd[1][i][0];
    }
  }
  /*
  for (j=0; j<10; j++) 
    printf("%f %f\n",prd[0][0][j],prd[1][0][j]);  */

}

fixviews()
{
  int i, j;
  float p0r, p0i,prr, pri, cc, dd, tmp;


  /* remove mag and phase variations (e.g. Hu's method) */
  p0r = p0i = 0.;
  for (i=0; i<npr; i++) {
    for (j=0; j<= samp1; j++) { 
      p0r += prd[0][i][j];
      p0i += prd[1][i][j];
    }
  }
  p0r /= (1.0 * npr);
  p0i /= (1.0 * npr);

  /* printf("Correction angle = "); */ 
  for (i=0; i<npr; i++) {
    prr = pri = 0.;
    for (j=0; j<= samp1; j++) { 
      prr += prd[0][i][j];
      pri += prd[1][i][j];
    }
    /*
    printf("%f %f\n",hypot(pri,prr),atan2(pri,prr));   */
    /* phase correction angle */
    tmp = atan2(p0i,p0r) - atan2(pri,prr); 
    if (tmp >  PI) tmp -= 2.0*PI;
    if (tmp < -PI) tmp += 2.0*PI;
    cc = cos(tmp*phasefact);
    dd = sin(tmp*phasefact);
    /* magnitude correction term 
    tmp = hypot(p0r,p0i) / (hypot(prr,pri));
    cc *= tmp*magfact + (1.0 - magfact);
    dd *= tmp*magfact + (1.0 - magfact);  */
    /* printf("%4.2f, ",atan2(dd,cc)*180.0/PI); */ 
    for (j=0; j< ndat; j++) { 
      prr = prd[0][i][j];
      pri = prd[1][i][j];
      tmp = cc*prr - dd*pri;
      prd[1][i][j] = dd*prr + cc*pri;
      prd[0][i][j] = tmp;
    }
  }
  /* printf("\n"); */
}

load_ft1()
{
  int i, j, lx, ly;
  float kx,ky,w,pr,pi,mkr,dkx,dky,dwin,w2,wx;
  float tmp1,tmpd,rot1,rot2,rotfact;
  int imnum,tmpslice,ie;
  int lxmn,lxmx,lymn,lymx;

  if (reg_flg) {
    ie = slnum+1;
    imnum = imoffset + rawnum*nphases + phnum; /* +1 if not Mark's numbering */
    if (sliceorder)
      tmpslice = ie;
    else
      tmpslice = (ie <= (nslices+1)/2) ? (ie*2 - 1) : ((ie-(nslices+1)/2)*2);
    tmpslice -= 1; /* for Mark's numbering scheme */
    if (rev_flg) tmpslice = nslices - tmpslice - 1;
    /* printf("%d %d %f %f %f\n",imnum,tmpslice,regxs[tmpslice][imnum],regys[tmpslice][imnum],regrot[tmpslice][imnum]); */
  }

  /* convolution function ---- kaiser-bessel */
  for (i=0; i < NWEIGHTS+1; i++)
    weight[i] = kaiser((float)NWEIGHTS,gridb,(float)(i-NWEIGHTS/2));
  tmp1 = (24.7/weight[NWEIGHTS/2]);
  for (i=0; i < NWEIGHTS+1; i++)
  {
    weight[i] *= tmp1;
    /* printf("%f\n",weight[i]); */
  }

  for (i=0; i<nim*OVER; i++)
    for (j=0; j<nim*OVER; j++) {
      sampim[i][j] = 0.0;
      grim[0][i][j] = 0.0;
      grim[1][i][j] = 0.0;
      ws[i][j] = 0.0;
  }

  /* dwin = 0.5*wind*OVER; */
  maxk = 0.;
  dwin = wind;
  rotfact = 2.0*PI/nim/OVER;
  for (i=0; i<npr; i++) {
    for (j=0; j< (ndat-samp1-1); j++) { 
      w2 = kdens[j];
      if (rsp_flg==0) {
        kx = t2k[0][i][j];
        ky = t2k[1][i][j];
      }
      else {
        kx = -t2k[0][i][j];
        ky = -t2k[1][i][j];
      }
      if ( (mkr = hypot(kx,ky)) > maxk) maxk = mkr;
      pr = prd[0][i][j+samp1];
      pi = prd[1][i][j+samp1];
      /* printf("%f %f %f %f %f\n",kx,ky,w2,pr,pi); DCN  */
      dkx = (factxx*kx + factxy*ky)*OVER; 
      dky = (factyx*kx + factyy*ky)*OVER; 
      if (reg_flg) { 
	rot1 = cos(rotfact*dkx*(pix_shifth*nim - regxs[tmpslice][imnum] + lrshift));
	rot2 = -sin(rotfact*dkx*(pix_shifth*nim - regxs[tmpslice][imnum] + lrshift));
	tmpd = pr*rot1 + pi*rot2;
	pi = pi*rot1 - pr*rot2;
	pr = tmpd;

	rot1 = cos(rotfact*dky*(pix_shiftv*nim - regys[tmpslice][imnum] + tbshift));
	rot2 = -sin(rotfact*dky*(pix_shiftv*nim - regys[tmpslice][imnum] + tbshift));
	tmpd = pr*rot1 + pi*rot2;
	pi = pi*rot1 - pr*rot2;
	pr = tmpd;

	rot1 = cos(-PI/180.0*regrot[tmpslice][imnum]);
	rot2 = -sin(-PI/180.0*regrot[tmpslice][imnum]);
	tmpd = dkx*rot1 + dky*rot2;
	dky = dky*rot1 - dkx*rot2;
	dkx = tmpd; 
      }
      else if (lrshift || tbshift || loc_flg) { 
	rot1 = cos(rotfact*dkx*(pix_shifth*nim + lrshift));
	rot2 = -sin(rotfact*dkx*(pix_shifth*nim + lrshift));
	tmpd = pr*rot1 + pi*rot2;
	pi = pi*rot1 - pr*rot2;
	pr = tmpd;

	rot1 = cos(rotfact*dky*(pix_shiftv*nim + tbshift));
	rot2 = -sin(rotfact*dky*(pix_shiftv*nim + tbshift));
	tmpd = pr*rot1 + pi*rot2;
	pi = pi*rot1 - pr*rot2;
	pr = tmpd;
      }
      dkx -= refl[2]*j*OVER; 
      dky -= refl[1]*j*OVER; 
      /* now do rotations and zooms */
      tmpd = dkx*factxxz + dky*factxyz;
      dky = factyxz*dkx + factyyz*dky; 
      dkx = tmpd; 
      dkx += nim/2*OVER; 
      dky += nim/2*OVER; 
      for (lx = ceil(dkx-dwin); lx<=floor(dkx+dwin); lx++) {
        if ((lx<0) || (lx>=OVER*nim)) continue;
	wx = weight[ (int) rint((((dkx-lx)/dwin)*NWEIGHTS/2 )) + NWEIGHTS/2] *w2;
        for (ly = ceil(dky-dwin); ly<=floor(dky+dwin); ly++) {
          if ((ly<0) || (ly>=OVER*nim)) continue;
	  w = wx*weight[ (int) rint((((dky-ly)/dwin)*NWEIGHTS/2 )) + NWEIGHTS/2];
	  sampim[lx][ly] += w * j;
	  grim[0][lx][ly] += pr * w;
	  grim[1][lx][ly] += pi * w;
	  ws[lx][ly] += w;
	}
      }
    }
  }
  
  for (i=0; i<nim*OVER; i++)
    for (j=0; j<nim*OVER; j++) {
      sampim[i][j] /= (0.001+ws[i][j]);
    }
}

load_ft2()
{
  int i, j, lx, ly;
  float kx,ky,w,pr,pi,mkr,dkx,dky,dwin,w2,wx;
  float tmpd,rot1,rot2,rotfact;
  int lxmn,lxmx,lymn,lymx;


  for (i=0; i<nim*OVER; i++)
    for (j=0; j<nim*OVER; j++) {
      grim[0][i][j] = 0.0;
      grim[1][i][j] = 0.0;
  }

  /* dwin = 0.5*wind*OVER; */
  dwin = wind;
  rotfact = 2.0*PI/nim/OVER;
  for (i=0; i<npr; i++) {
    for (j=0; j< (ndat-samp1-1); j++) { 
      w2 = kdens[j];
      if (rsp_flg==0) {
        kx = t2k[0][i][j];
        ky = t2k[1][i][j];
      }
      else {
        kx = -t2k[0][i][j];
        ky = -t2k[1][i][j];
      }
      pr = prd[0][i][j+samp1];
      pi = prd[1][i][j+samp1];
      dkx = (factxx*kx + factxy*ky)*OVER; 
      dky = (factyx*kx + factyy*ky)*OVER; 
      if (lrshift || tbshift || loc_flg) { 
	rot1 = cos(rotfact*dkx*(pix_shifth*nim + lrshift));
	rot2 = -sin(rotfact*dkx*(pix_shifth*nim + lrshift));
	tmpd = pr*rot1 + pi*rot2;
	pi = pi*rot1 - pr*rot2;
	pr = tmpd;

	rot1 = cos(rotfact*dky*(pix_shiftv*nim + tbshift));
	rot2 = -sin(rotfact*dky*(pix_shiftv*nim + tbshift));
	tmpd = pr*rot1 + pi*rot2;
	pi = pi*rot1 - pr*rot2;
	pr = tmpd;
      }
      /* now do rotations and zooms */
      tmpd = dkx*factxxz + dky*factxyz;
      dky = factyxz*dkx + factyyz*dky; 
      dkx = tmpd; 
      dkx += nim/2*OVER; 
      dky += nim/2*OVER; 
      for (lx = ceil(dkx-dwin); lx<=floor(dkx+dwin); lx++) {
        if ((lx<0) || (lx>=OVER*nim)) continue;
	wx = weight[ (int) rint((((dkx-lx)/dwin)*NWEIGHTS/2 )) + NWEIGHTS/2] *w2;
        for (ly = ceil(dky-dwin); ly<=floor(dky+dwin); ly++) {
          if ((ly<0) || (ly>=OVER*nim)) continue;
	  w = wx*weight[ (int) rint((((dky-ly)/dwin)*NWEIGHTS/2 )) + NWEIGHTS/2];
	  grim[0][lx][ly] += pr * w;
	  grim[1][lx][ly] += pi * w;
	}
      }
    }
  }
}

fermi_filt1()
{
  float ww,rad;
  int i,j,offs;
  rad = nim/2*OVER;
  if ( maxk*OVER/zoomer < rad )  rad = maxk*OVER/zoomer;
  /* printf("rad = %f ",rad); */
  offs = nim/2*OVER;
  for(i=0; i<nim*OVER; i++) 
    for (j=0; j<nim*OVER; j++) {
      ws[i][j] = 1. / (1. + exp((hypot((float)(offs-i),(float)(offs-j)) - rad) *6.4 /rad));
/*    ws[i][j] = 1. / ((1. + exp((hypot((float)(offs-i),(float)(offs-j)) - rad) *6.4 /rad))*(0.001+ws[i][j]));  - with density post-comp */
    }
}

fermi_filt2()
{
  int i,j;
  for(i=0; i<nim*OVER; i++) 
    for (j=0; j<nim*OVER; j++) {
      grim[0][i][j] *= ws[i][j];
      grim[1][i][j] *= ws[i][j];
    }
}


copy_dat()
{
  int i,j;

  for(i=0; i<nim*OVER; i++) 
    for (j=0; j<nim*OVER; j++) {
      im[0][i][j] = grim[0][i][j];
      im[1][i][j] = grim[1][i][j];
    }
}

wind_dat(midsamp,sampwind)
int midsamp,sampwind;
{
  float ww,rad;
  float stndat, enndat;
  int i,j,offs;

  stndat = midsamp - sampwind;
  enndat = midsamp + sampwind;

  for(i=0; i<nim*OVER; i++) 
    for (j=0; j<nim*OVER; j++) {
      rad = sampim[i][j];
      if ((rad > stndat) && (rad < enndat))
        ww = .5*(1.0 + cos(PI*(rad-midsamp)/sampwind));
      else ww = 0.;
      im[0][i][j] = ww*grim[0][i][j];
      im[1][i][j] = ww*grim[1][i][j];
    }
}

flipx_image()
{
  int i,j, n;

  n = nim*OVER; 
  for(i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      w1[0][j] = im[0][j][i];
      w1[1][j] = im[1][j][i];
    }
    for (j=0; j<n; j++) {
      im[0][j][i] = w1[0][n-j-1];
      im[1][j][i] = w1[1][n-j-1];
    }
  }
}

flipy_image()
{
  int i,j, n;

  n = nim*OVER; 
  for(i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      w1[0][j] = im[0][i][j];
      w1[1][j] = im[1][i][j];
    }
    for (j=0; j<n; j++) {
      im[0][i][j] = w1[0][n-j-1];
      im[1][i][j] = w1[1][n-j-1];
    }
  }
}

ft_image()
{
  int i,j, fftd, n,offs,offs2;

  n = nim*OVER; fftd = 1;
  offs = n/2;
  offs2 = (OVER-1)*nim/2;
  for(i=0; i<nim*OVER; i++) {
    for (j=0; j<offs; j++) {
      w1[0][j] = im[0][i][j+offs];
      w1[1][j] = im[1][i][j+offs];
      w1[0][j+offs] = im[0][i][j];
      w1[1][j+offs] = im[1][i][j];
    }
    fft(&fftd, &n, w1[0], w1[1]);
    for (j=0; j<offs; j++) {
      im[0][i][j+offs] = w1[0][j];
      im[1][i][j+offs] = w1[1][j];
      im[0][i][j] = w1[0][j+offs];
      im[1][i][j] = w1[1][j+offs];
    }
  }
  for(j=0; j<nim; j++) {
    for (i=0; i<offs; i++) {
      w1[0][i] = im[0][i+offs][j+offs2];
      w1[1][i] = im[1][i+offs][j+offs2];
      w1[0][i+offs] = im[0][i][j+offs2];
      w1[1][i+offs] = im[1][i][j+offs2];
    }
    fft(&fftd, &n, w1[0], w1[1]);
    for (i=0; i<offs; i++) {
      im[0][i+offs][j+offs2] = w1[0][i];
      im[1][i+offs][j+offs2] = w1[1][i];
      im[0][i][j+offs2] = w1[0][i+offs];
      im[1][i][j+offs2] = w1[1][i+offs];
    }
  }
}

make_final(segnum)
int segnum;
{
  int i, j,offs;
  float phfact;

  if (rsp_flg==0) 
    phfact = -segnum;
  else
    phfact = segnum;
  offs = (OVER-1)*nim/2;
  if (segnum == 0)
    for (i=0; i<nim; i++)
      for (j=0; j<nim; j++)
      {
        finalim[0][i][j] = im[0][i+offs][j+offs];
        finalim[1][i][j] = im[1][i+offs][j+offs];
      }
  else
    for (i=0; i<nim; i++)
      for (j=0; j<nim; j++)
      {
        finalim[0][i][j] += 
	      im[0][i+offs][j+offs]*cos(phfact*refim[i][j]*SEGSIZE/mapdel) -
	      im[1][i+offs][j+offs]*sin(phfact*refim[i][j]*SEGSIZE/mapdel);
        finalim[1][i][j] += 
	      im[0][i+offs][j+offs]*sin(phfact*refim[i][j]*SEGSIZE/mapdel) +
	      im[1][i+offs][j+offs]*cos(phfact*refim[i][j]*SEGSIZE/mapdel);
      }
  if (segnum == nsegs)
    for (i=0; i<nim; i++)
      for (j=0; j<nim; j++)
      {
	      im[0][i+offs][j+offs] = finalim[0][i][j];
	      im[1][i+offs][j+offs] = finalim[1][i][j];
      } 
}

write_raw()
{
  int i, j, mx,offs,tmpslice,imnum;
  float imx, imi, imr, imm, wctmp, wcmax;
  short int bm[OVER*MXNIM];
  char fn[80];

  sprintf(fn,"raw.%d.%d.i",slnum+1,phnum+1);

  if (!(fo1=fopen(fn,"w"))) {
    fprintf(stderr,"Can't open %s!\n",fn); return(0); }

  imx = 0.0;
  for (i=0; i<OVER*nim; i++) {
    for (j=0; j<OVER*nim; j++) {
        imr = grim[0][j][i];
        imi = grim[1][j][i];
        imm = sqrt(imi*imi+imr*imr);
	if (imm > imx) imx = imm;
        bm[j] = imm;
    }
    fwrite(bm,OVER*nim,2,fo1);    
  }
  mx = imx;
  if (out_flg) printf("max=%3d",mx);
  fclose(fo1);
}

write_imagemc(nc,ie,nph,nraw)
int nc,ie,nph,nraw;
{
  int i, j, k, mx,mm,offs,tmpslice,imnum;
  float imx, imi, imr, imm, imt, wctmp, wcmax;
  float wc[MXNIM];
  short imtmp[MXNIM*MXNIM];
  short int bm[MXNIM];
  char fn[80];


  imnum = imoffset + nraw*nphases + nph;
  if (sliceorder)
    tmpslice = ie;
  else
    tmpslice = (ie <= (nslices+1)/2) ? (ie*2 - 1) : ((ie-(nslices+1)/2)*2);
  if (rev_flg) tmpslice = nslices - tmpslice + 1;

  sprintf(fn,"sl%1d.%03d",tmpslice,imnum);
  
  if ((nc != 1)||(L2_flg)) {
    if (!(fo1=fopen(fn,"r"))) {
      fprintf(stderr,"Can't open %s!\n",fn); return(0); }

    if (fread(imtmp,sizeof(short),nim*nim,fo1) != nim*nim ) {
      fprintf(stderr,"Error reading %s!\n",fn); return(0); }

    fclose(fo1);
  } else {
    for (i=0; i<nim*nim; i++) 
      imtmp[i] = 0;
  }

  if (!(fo1=fopen(fn,"w"))) {
    fprintf(stderr,"Can't open %s!\n",fn); return(0); }

  /* kaiser-bessel correction */
  wcmax = sinh(sqrt(gridb*gridb))/sqrt(gridb*gridb);

  for (i=0; i<nim; i++)
  {
    wctmp = PI*PI*gridl*gridl*(i-nim/2)*(i-nim/2)/(nim*nim*OVER*OVER/4) - gridb*gridb;
    if(wctmp == 0.)
      wc[i] = 1.0*wcmax;
    else if(wctmp < 0.)
      wc[i] = sqrt(-wctmp)/sinh(sqrt(-wctmp))*wcmax;
    else
      wc[i] = sqrt(wctmp)/sin(sqrt(wctmp))*wcmax;
  }

  imx = 0.0;
  offs = (OVER-1)*nim/2;
  for (i=0; i<nim; i++) {
    for (j=0; j<nim; j++) {
      if (hypot((float)(i-nim/2),(float)(j-nim/2)) < .51*nim)
      {
	imt = imtmp[i*nim + j];
	imt *= (nim*nim*OVER*OVER)/outscale;
        imr = im[0][-j+nim+offs-1][-i+nim+offs-1]*wc[i]*wc[j];
        imi = im[1][-j+nim+offs-1][-i+nim+offs-1]*wc[i]*wc[j];
        imm = sqrt(imi*imi+imr*imr + imt*imt);
	if (imm > imx) imx = imm;
        bm[j] = outscale*imm/(nim*nim*OVER*OVER);
      }
      else
        bm[j] = 0;
    }
    fwrite(bm,nim,2,fo1);    
  }
  mx = imx*outscale/(nim*nim*OVER*OVER);
  if (out_flg) printf("max=%3d",mx);
  fclose(fo1);
}

write_anal(ie)
int ie;
{
  int i, j, mx,offs,tmpslice,imnum;
  float imx, imi, imr, imm, wctmp, wcmax;
  float wc[MXNIM];
  double outnum;

  /* kaiser-bessel correction */
  wcmax = sinh(sqrt(gridb*gridb))/sqrt(gridb*gridb);

  for (i=0; i<nim; i++)
  {
    wctmp = PI*PI*gridl*gridl*(i-nim/2)*(i-nim/2)/(nim*nim*OVER*OVER/4) - gridb*gridb;
    if(wctmp == 0.)
      wc[i] = 1.0*wcmax;
    else if(wctmp < 0.)
      wc[i] = sqrt(-wctmp)/sinh(sqrt(-wctmp))*wcmax;
    else
      wc[i] = sqrt(wctmp)/sin(sqrt(wctmp))*wcmax;
  }

  if (sliceorder)
    tmpslice = ie;
  else
    tmpslice = (ie <= (nslices+1)/2) ? (ie*2 - 1) : ((ie-(nslices+1)/2)*2);
  if (rev_flg) tmpslice = nslices - tmpslice + 1;
  tmpslice--;

  offs = (OVER-1)*nim/2;
  for (i=0; i<nim; i++) {
    for (j=0; j<nim; j++) {
      if (hypot((float)(i-nim/2),(float)(j-nim/2)) < .51*nim)
      {
       imr = im[0][-j+nim+offs-1][-i+nim+offs-1]*wc[i]*wc[j];
       imi = im[1][-j+nim+offs-1][-i+nim+offs-1]*wc[i]*wc[j];
       outnum = outscale*sqrt(imi*imi+imr*imr)/(nim*nim*OVER*OVER);
       if (isnan(outnum))
         outbuf[tmpslice][i][j] = 0;
       else
         outbuf[tmpslice][i][j] = outnum;
      }
      else
        outbuf[tmpslice][i][j] = 0;
    }
  }
}

write_anal_all(nc,nph,nraw)
int nc,nph,nraw;
{
  int i, j, k, mx,mn,offs,tmpslice,imnum;
  char fn[80];
  struct hdr hdr;
  float imr,imi;

  imnum = imoffset + nraw*nphases + nph;

  /* output individual coils */
  if (allcoils_flg) {
    sprintf(fn,"vol_e%d_%s_%04d_%d.img",exam,datestr,imnum,nc); 
    if (!(fo1=fopen(fn,"w"))) {
      fprintf(stderr,"Can't open %s!\n",fn); return(0); }
    for (i=0; i<nslices; i++) 
      for (j=0; j<nim; j++) 
        fwrite(outbuf[i][j],nim,2,fo1);    
    fclose(fo1);
  }

  sprintf(fn,"vol_e%d_%s_%04d.img",exam,datestr,imnum); 
  /* sprintf(fn,"vol_e%d_%04d.img",exam,imnum); */

  /* section to do mc and spiral in/out combintation */
  if ((nc != 1)||(L2_flg)) {
    if (!(fo1=fopen(fn,"r"))) {
      fprintf(stderr,"Can't open %s!\n",fn); return(0); }

    for (i=0; i<nslices; i++) 
      for (j=0; j<nim; j++) {
        if (fread(outbuf2[i][j],sizeof(short),nim,fo1) != nim ) {
          fprintf(stderr,"Error reading %s (%d,%d) !\n",fn,i,j); return(0); }
        for (k=0; k<nim; k++) {
          imr = outbuf[i][j][k];
          imi = outbuf2[i][j][k];
	  outbuf[i][j][k] = sqrt(imr*imr + imi*imi);
        }
      }
    fclose(fo1);
  } 

  if (!(fo1=fopen(fn,"w"))) {
    fprintf(stderr,"Can't open %s!\n",fn); return(0); }

  mn = mx = outbuf[0][0][0];
  for (i=0; i<nslices; i++) {
    for (j=0; j<nim; j++) {
      for (k=0; k<nim; k++) {
	if (outbuf[i][j][k] > mx) mx = outbuf[i][j][k];
	if (outbuf[i][j][k] < mn) mn = outbuf[i][j][k];
      }
      fwrite(outbuf[i][j],nim,2,fo1);    
    }
  }
  if (out_flg) printf("min/max=%d/%d\n",mn,mx);
  fclose(fo1);

  /* now do analyze header */
  sprintf(fn,"vol_e%d_%s_%04d.hdr",exam,datestr,imnum); 
  /* sprintf(fn,"vol_e%d_%04d.hdr",exam,imnum); */

  if (!(fo1=fopen(fn,"w"))) {
    fprintf(stderr,"Can't open %s!\n",fn); return(0); }

  {
      /*Initialize hdr to zeros to allow compatibility with newest*/
      /*versions of ANALYZE*/
      /*With thanks to Joel T. Lee, Psychiatry PET, VA, Minneapolis*/
      char *ptr,*endptr;
      for(ptr=(char *)&hdr,endptr=ptr+sizeof(hdr);ptr<endptr;ptr++)*ptr=0;
  }
	/*Copy required data into header struct*/
	hdr.bits=16; 
	hdr.datatype=4;		/*short int*/
	hdr.dims=4;
	hdr.x_dim=nim;
	hdr.y_dim=nim;
	hdr.z_dim=nslices;
	hdr.t_dim=1;
	hdr.x_size=10.0*opfov/nim;
	hdr.y_size=10.0*opfov/nim;
	hdr.z_size=slthick;
	hdr.glmax=mx;			/*global max*/
	hdr.glmin=mn;			/*global min*/
	hdr.sizeof_hdr=sizeof(struct hdr);	/*standard=348*/
	hdr.extents=16384;
	hdr.regular='r';

	/*Write out header*/
	if(fwrite((char *)&hdr,sizeof(struct hdr),1,fo1)!=1){
		printf("file write error writing header %s\n",fn);
		fclose(fo1);
		return 0;
	}
	fclose(fo1);

	if (allcoils_flg) {
          sprintf(fn,"vol_e%d_%s_%04d_%d.hdr",exam,datestr,imnum,nc); 
          if (!(fo1=fopen(fn,"w"))) {
            fprintf(stderr,"Can't open %s!\n",fn); return(0); }
	  if(fwrite((char *)&hdr,sizeof(struct hdr),1,fo1)!=1){
		printf("file write error writing header %s\n",fn);
		fclose(fo1);
		return 0;
	  }
	  fclose(fo1);
	}
}

write_image(nc,ie,nph,nraw)
int nc,ie,nph,nraw;
{
  int i, j, mx,offs,tmpslice,imnum;
  float imx, imi, imr, imm, wctmp, wcmax;
  float wc[MXNIM];
  short int bm[MXNIM];
  char fn[80];

  /* kaiser-bessel correction */
  wcmax = sinh(sqrt(gridb*gridb))/sqrt(gridb*gridb);

  for (i=0; i<nim; i++)
  {
    wctmp = PI*PI*gridl*gridl*(i-nim/2)*(i-nim/2)/(nim*nim*OVER*OVER/4) - gridb*gridb;
    if(wctmp == 0.)
      wc[i] = 1.0*wcmax;
    else if(wctmp < 0.)
      wc[i] = sqrt(-wctmp)/sinh(sqrt(-wctmp))*wcmax;
    else
      wc[i] = sqrt(wctmp)/sin(sqrt(wctmp))*wcmax;
  }

  imnum = imoffset + nraw*nphases + nph;
  if (sliceorder)
    tmpslice = ie;
  else
    tmpslice = (ie <= (nslices+1)/2) ? (ie*2 - 1) : ((ie-(nslices+1)/2)*2);
  if (rev_flg) tmpslice = nslices - tmpslice + 1;

  sprintf(fn,"sl%1d.%1d.%03d",tmpslice,nc,imnum);

  if (!(fo1=fopen(fn,"w"))) {
    fprintf(stderr,"Can't open %s!\n",fn); return(0); }

  imx = 0.0;
  offs = (OVER-1)*nim/2;
  for (i=0; i<nim; i++) {
    for (j=0; j<nim; j++) {
      if (hypot((float)(i-nim/2),(float)(j-nim/2)) < .51*nim)
      {
        imr = im[0][-j+nim+offs-1][-i+nim+offs-1]*wc[i]*wc[j];
        imi = im[1][-j+nim+offs-1][-i+nim+offs-1]*wc[i]*wc[j];
        imm = sqrt(imi*imi+imr*imr);
	if (imm > imx) imx = imm;
        bm[j] = outscale*imm/(nim*nim*OVER*OVER);
      }
      else
        bm[j] = 0;
    }
    fwrite(bm,nim,2,fo1);    
  }
  mx = imx*outscale/(nim*nim*OVER*OVER);
  if (out_flg) printf("max=%3d",mx);
  fclose(fo1);
}

write_image_phs(ie,nph,nraw)
int ie,nph,nraw;
{
  int i,j,offs,tmpslice,imnum;
  float imi, imr;
  short int bm[MXNIM];
  char fn[80];

  imnum = imoffset + nraw*nphases + nph;
  if (sliceorder)
    tmpslice = ie;
  else
    tmpslice = (ie <= (nslices+1)/2) ? (ie*2 - 1) : ((ie-(nslices+1)/2)*2);
  if (rev_flg) tmpslice = nslices - tmpslice + 1;

  sprintf(fn,"sl%1d.phs.%03d",tmpslice,imnum);

  if (!(fo1=fopen(fn,"w"))) {
    fprintf(stderr,"Can't open %s!\n",fn); return(0); }

  offs = (OVER-1)*nim/2;
  for (i=0; i<nim; i++) {
    for (j=0; j<nim; j++) {
      if (hypot((float)(i-nim/2),(float)(j-nim/2)) < .51*nim)
      {
        imr = im[0][-j+nim+offs-1][-i+nim+offs-1];
        imi = im[1][-j+nim+offs-1][-i+nim+offs-1];
        bm[j] = 1000 * atan2(imi,imr);
      }
      else
        bm[j] = 0;
    }
    fwrite(bm,nim,2,fo1);    
  }
  fclose(fo1);
}

make_phase1()
{
  int i, j,offs;

  offs = (OVER-1)*nim/2;
  for (i=0; i<nim; i++)
    for (j=0; j<nim; j++)
      refim[i][j] = atan2(im[1][i+offs][j+offs],im[0][i+offs][j+offs]);
}

make_phase2()
{
  int i, j,offs;
  float phaseoff;

  phaseoff = -2.0*PI*(mapdel/samptime)*pt/ndat; 
  offs = (OVER-1)*nim/2;
  for (i=0; i<nim; i++)
    for (j=0; j<nim; j++)
    {
      refimmag[i][j] = hypot(im[1][i+offs][j+offs],im[0][i+offs][j+offs]);
      refim[i][j] -= (atan2(im[1][i+offs][j+offs],im[0][i+offs][j+offs]) 
			- phaseoff);
    }
}

make_phase1_filt()
{
  int i, j, k, offs;
  int hfsize;
  float fbuf[2][OVER*MXNIM];
  double sumr,sumi;

  hfsize = (fsize-1)/2;
  offs = (OVER-1)*nim/2;

  for (i=0; i<nim; i++)
  {
    for (j= MAX(0,(offs-hfsize)); j<= MIN((nim+offs+hfsize),OVER*MXNIM-1); j++)
    {
      fbuf[0][j] = 0.;
      fbuf[1][j] = 0.;
      for (k= MAX(0,(i+offs-hfsize)); k<= MIN((i+offs+hfsize),OVER*MXNIM-1); k++)
      {
        fbuf[0][j] += im[0][k][j];
        fbuf[1][j] += im[1][k][j];
      }
    }
    for (j=0; j<nim; j++)
    {
      sumr = sumi = 0.;
      for (k= MAX(0,(j+offs-hfsize)); k<= MIN((j+offs+hfsize),OVER*MXNIM-1); k++)
      {
        sumr += fbuf[0][k];
        sumi += fbuf[1][k];
      }
      refim[i][j] = atan2(sumi,sumr);
    }
  }
}

make_phase2_filt()
{
  int i, j, k, offs;
  int hfsize;
  float fbuf[2][OVER*MXNIM];
  double sumr,sumi,phaseoff;

  /* phaseoff is value to offset freq map if phase twist (pt) != 0 */
  phaseoff = -2.0*PI*(mapdel/samptime)*pt/ndat; 
  hfsize = (fsize-1)/2;
  offs = (OVER-1)*nim/2;

  for (i=0; i<nim; i++)
  {
    for (j= MAX(0,(offs-hfsize)); j<= MIN((nim+offs+hfsize),OVER*MXNIM-1); j++)
    {
      fbuf[0][j] = 0.;
      fbuf[1][j] = 0.;
      for (k= MAX(0,(i+offs-hfsize)); k<= MIN((i+offs+hfsize),OVER*MXNIM-1); k++)
      {
        fbuf[0][j] += im[0][k][j];
        fbuf[1][j] += im[1][k][j];
      }
    }
    for (j=0; j<nim; j++)
    {
      sumr = sumi = 0.;
      for (k= MAX(0,(j+offs-hfsize)); k<= MIN((j+offs+hfsize),OVER*MXNIM-1); k++)
      {
        sumr += fbuf[0][k];
        sumi += fbuf[1][k];
      }
      refimmag[i][j] = hypot(sumi,sumr);
      /* refim[i][j] = atan2(sumi,sumr); */
      refim[i][j] -= (atan2(sumi,sumr) - phaseoff); 
    }
  }
}

write_ref(ie)
int ie;
{
  int i, j;
  float bm[MXNIM];
  char fn[80];
  float max;
  double tmpxr, tmpyr, tmpxi, tmpyi, tmpc, tmpcr, tmpci;
  float out[3];

  if (linmap_flg) {

    max = 0.0;
    for (i=0; i<nim; i++) 
      for (j=0; j<nim; j++) 
        if (hypot((float)(i-nim/2),(float)(j-nim/2)) < nim/2) 
	  if (refimmag[i][j] > max) max = refimmag[i][j];
    max /= 4.0;

    tmpxr = tmpyr = tmpxi = tmpyi = tmpcr = tmpci = 0.0;
    /* find linear terms */
    for (i=0; i<nim-1; i++) 
      for (j=0; j<nim-1; j++) 
	if (refimmag[i][j] > max) {
	  tmpxr += refimmag[i][j]*refimmag[i][j+1]
		   *cos(refim[i][j]-refim[i][j+1]);
	  tmpxi += refimmag[i][j]*refimmag[i][j+1]
		   *sin(refim[i][j]-refim[i][j+1]);
	  tmpyr += refimmag[i][j]*refimmag[i+1][j]
		   *cos(refim[i][j]-refim[i+1][j]);
	  tmpyi += refimmag[i][j]*refimmag[i+1][j]
		   *sin(refim[i][j]-refim[i+1][j]);
        }

    tmpxr = -atan2(tmpxi,tmpxr);
    tmpyr = -atan2(tmpyi,tmpyr);
    /* uncomment to zero linear terms */
    /* tmpxr = tmpyr = 0.; */ 

    /* find constant term with linears removed */
    for (i=0; i<nim-1; i++) 
      for (j=0; j<nim-1; j++) 
	if (refimmag[i][j] > max) {
	  tmpcr += refimmag[i][j]
	      *cos(refim[i][j] - tmpxr*(j-nim/2) - tmpyr*(i-nim/2));
	  tmpci += refimmag[i][j]
	      *sin(refim[i][j] - tmpxr*(j-nim/2) - tmpyr*(i-nim/2));
        }

    tmpc = atan2(tmpci,tmpcr);
    /* 
    out[0] = tmpc;
    out[1] = tmpxr*nim;
    out[2] = tmpyr*nim;
    printf("image = %d, constant = %f, x-grad = %f, y-grad = %f\n",
      ie,out[0],out[1],out[2]); */

    out[0] = tmpc*samptime/mapdel;
    out[1] = tmpxr*nim*samptime/mapdel/2.0/PI;
    out[2] = tmpyr*nim*samptime/mapdel/2.0/PI;

    sprintf(fn,"refl.s%1d",ie);
    if (!(fo1=fopen(fn,"w"))) {
      fprintf(stderr,"Can't open %s!\n",fn); return(0); }
    fwrite(out,3,sizeof(*out),fo1);    
    fclose(fo1);

    for (i=0; i<nim-1; i++) 
      for (j=0; j<nim-1; j++) {
	  tmpcr = refimmag[i][j]
	      *cos(refim[i][j] - tmpxr*(j-nim/2) - tmpyr*(i-nim/2) - tmpc);
	  tmpci = refimmag[i][j]
	      *sin(refim[i][j] - tmpxr*(j-nim/2) - tmpyr*(i-nim/2) - tmpc);
          refim[i][j]= atan2(tmpci,tmpcr);
        }
  }

  if (map_flg) {
    sprintf(fn,"ref.%d.s%1d",nim,ie);
    if (!(fo1=fopen(fn,"w"))) {
      fprintf(stderr,"Can't open %s!\n",fn); return(0); }

    for (i=0; i<nim; i++) {
      for (j=0; j<nim; j++) {
        if (refim[i][j] > PI) refim[i][j] -= 2*PI;
        if (refim[i][j] < -PI) refim[i][j] += 2*PI;
        if (refim[i][j] > 0.5*PI) refim[i][j] = 0.5*PI;
        if (refim[i][j] < -0.5*PI) refim[i][j] = -0.5*PI;
        bm[j] = refim[i][j];
      }
      fwrite(bm,nim,sizeof(*bm),fo1);    
    }
    fclose(fo1);
  }
}

write_ref_mc(ie,nc)
int ie,nc;
{
  int i, j;
  float bm[MXNIM];
  char fn[80];

  sprintf(fn,"ref.%d.%d.s%1d",nim,nc,ie);
  if (!(fo1=fopen(fn,"w"))) {
    fprintf(stderr,"Can't open %s!\n",fn); return(0); }

  for (i=0; i<nim; i++) {
    for (j=0; j<nim; j++) {
        if (refim[i][j] > PI) refim[i][j] -= 2*PI;
        if (refim[i][j] < -PI) refim[i][j] += 2*PI;
        bm[j] = refim[i][j];
    }
    fwrite(bm,nim,sizeof(*bm),fo1);    
  }
  fclose(fo1);

  sprintf(fn,"refm.%d.%d.s%1d",nim,nc,ie);
  if (!(fo1=fopen(fn,"w"))) {
    fprintf(stderr,"Can't open %s!\n",fn); return(0); }

  for (i=0; i<nim; i++) {
    for (j=0; j<nim; j++) {
        bm[j] = refimmag[i][j];
    }
    fwrite(bm,nim,sizeof(*bm),fo1);    
  }
  fclose(fo1);
}

write_ref_mcf(ie)
int ie;
{
  int i, j, k;
  float bm[MXNIM];
  char fn[80];
  float maptmp[MXNIM][MXNIM];

  for (i=0; i<nim; i++) {
    for (j=0; j<nim; j++) {
	refim[i][j] = 0;
	refimmag[i][j] = 0;
    }
  }

  for (k=0; k<ncoils; k++) {
    sprintf(fn,"ref.%d.%d.s%1d",nim,k,ie);
    if (!(fo1=fopen(fn,"r"))) {
      fprintf(stderr,"Can't open reference file %s!\n",fn); return(0); }

    for (i=0; i<nim; i++) {
      fread(bm,nim,sizeof(*bm),fo1);    
      for (j=0; j<nim; j++) {
	maptmp[i][j] = bm[j];
      }
    }
    fclose(fo1);
    unlink(fn);

    sprintf(fn,"refm.%d.%d.s%1d",nim,k,ie);
    if (!(fo1=fopen(fn,"r"))) {
      fprintf(stderr,"Can't open reference file %s!\n",fn); return(0); }

    for (i=0; i<nim; i++) {
      fread(bm,nim,sizeof(*bm),fo1);    
      for (j=0; j<nim; j++) {
	refim[i][j] += maptmp[i][j]*bm[j];
	refimmag[i][j] += bm[j];
      }
    }
    fclose(fo1);
    unlink(fn);
  }
  for (i=0; i<nim; i++)
    for (j=0; j<nim; j++) {
      refim[i][j] /= (refimmag[i][j] + 0.001);
    }
}

load_ref_all()
{
  int i, j,ie;
  float bm[MXNIM];
  char fn[80];

  for (ie = 0; ie < nslices; ie++) {
  if (lincor_flg) { 
    sprintf(fn,"refl.s%1d",ie);
    if (!(fo1=fopen(fn,"r"))) {
      fprintf(stderr,"Can't open reference file %s!\n",fn); return(0); }
    fread(reflbuf[ie],3,sizeof(*refl),fo1);    
    fclose(fo1);
    /* printf("image = %d, constant = %f, x-grad = %f, y-grad = %f\n",
      ie,refl[0],refl[1],refl[2]);  */
  }

  if (hc_flg) { 
    sprintf(fn,"ref.%d.s%1d",nim,ie);
    /* printf("%s ",fn); */
    if (!(fo1=fopen(fn,"r"))) {
      fprintf(stderr,"Can't open reference file %s!\n",fn); return(0); }

    for (i=0; i<nim; i++) {
      fread(bm,nim,sizeof(*bm),fo1);    
      for (j=0; j<nim; j++) {
	refbuf[ie][i][j] = bm[j];
      }
    }
    fclose(fo1);
  }
  }
}

load_ref(ie)
int ie;
{
  int i, j;

  if (lincor_flg) 
    for (i=0; i<nim; i++) 
      refl[i] = reflbuf[ie][i];

  if (hc_flg) 
    for (i=0; i<nim; i++) 
      for (j=0; j<nim; j++) 
	refim[i][j] = refbuf[ie][i][j];
}

float kaiser(l, b, u)
float l, b, u;
{
    float f;
    float bessi0();

    if (fabs(u) <= 0.5*fabs(l))
    {
        f = bessi0(b*sqrt(1.0-(2.0*u/l)*(2.0*u/l)));
        return f;
    }
    else
        return 0.0;
}

float bessi0(x)
float x;
{
	float ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}

void fft(fwd,n,xr, xi)
int *fwd, *n;
float xr[],xi[];
{
  void four1();
  four1(xr,xi, *n, *fwd);
}


#include <math.h>

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(rdat,idat,nn,isign)
float rdat[],idat[];
int nn,isign;
{
	int n,mmax,m,j,j1,istep,i,mmaxby2;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=0;i<nn;i++) {
		j1 = (j-1)/2;
		if (j1 > i) {
			SWAP(rdat[j1],rdat[i]);
			SWAP(idat[j1],idat[i]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=2*mmax;
		mmaxby2 = mmax/2;
		theta=6.28318530717959/(isign*mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=0;m<mmaxby2;m++) {
			for (i=m;i<nn;i+=mmax) {
				j=i+mmaxby2;
				tempr=wr*rdat[j]-wi*idat[j];
				tempi=wr*idat[j]+wi*rdat[j];
				rdat[j]=rdat[i]-tempr;
				idat[j]=idat[i]-tempi;
				rdat[i] += tempr;
				idat[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

#undef SWAP


/* routine to get floats out of GE header - needed because
   header is not word aligned */
float headflt(addr)
char *addr;
{
   float fltval;
   memcpy((char*)&fltval,addr,4);
   return(fltval);
}

int headint(addr)
char *addr;
{
   int intval;
   memcpy((char*)&intval,addr,4);
   return(intval);
}

short headshrt(addr)
char *addr;
{
   short shrtval;
   memcpy((char*)&shrtval,addr,4);
   return(shrtval);
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void loc_calc(p1,p2,p3,rot,trans)
float *p1,*p2,*p3;
int rot,trans;
{
  /* 
    top left     - p1
    bottom left  - p2
    top right    - p3
    bottom right - p4
  */
  float p4[3],tbvec[3],lrvec[4],tempr;
  int i,n;

  /* calc missing corner (bottom right) */
  for (i=0; i<3; i++)
    p4[i] = p3[i] + (p2[i] - p1[i]);

  if (trans)
    for (i=0; i<3; i++) {
      SWAP(p2[i],p3[i]);
    }

  for (n=0; n<rot; n++)
    for (i=0; i<3; i++) {
      SWAP(p1[i],p3[i]);
      SWAP(p3[i],p4[i]);
      SWAP(p4[i],p2[i]);
    }

  /* reformat to match image header */
  for (i=0; i<3; i++) {
    SWAP(p2[i],p3[i]);
    SWAP(p3[i],p4[i]);
  }

}
#undef SWAP

