/********************** io_scidac.c **********************************/
/* MILC/QIO interface.  Produces an ILDG archivable file */
/* MIMD version 7 */
/* CD 11/2004
*/

#include "generic_includes.h"
#include <qio.h>
#include "../include/io_lat.h"
#include "../include/io_scidac.h"
#include <string.h>
#define LATDIM 4
void QIO_set_trelease(double t_in, double t_out);

/* Map QIO layout functions to MILC functions */

#ifndef QIO_HAS_EXTENDED_LAYOUT
#error Requires a QIO version that supports the extended layout for large files.
#endif

int qio_node_number(const int x[], void *arg){
  return node_number(x[0],x[1],x[2],x[3]);
}

QIO_Index qio_node_index(const int x[], void *arg){
  return node_index(x[0],x[1],x[2],x[3]);
}

static void qio_get_coords(int x[], int node, QIO_Index index, void *arg){
  /* For this node we have a table */
  if(node == this_node){
    x[0] = lattice[index].x;
    x[1] = lattice[index].y;
    x[2] = lattice[index].z;
    x[3] = lattice[index].t;
  }
  /* For other nodes we require the layout function */
  else
    get_coords( x, node, index );
}

QIO_Index qio_num_sites(int node, void *arg){
  return num_sites(node);
}

void build_qio_layout(QIO_Layout *layout){
  static int lattice_size[LATDIM];

  lattice_size[0] = nx;
  lattice_size[1] = ny;
  lattice_size[2] = nz;
  lattice_size[3] = nt;

  /* Members for smaller lattices are not used here. */
  layout->node_number = NULL;;
  layout->node_index  = NULL;
  layout->get_coords  = NULL;
  layout->num_sites   = 0;
  /* Members using the extended layout for large-lattice support */
  layout->node_number_ext = qio_node_number;
  layout->node_index_ext  = qio_node_index;
  layout->get_coords_ext  = qio_get_coords;
  layout->num_sites_ext   = qio_num_sites;
  layout->arg             = NULL;
  /* End of extended-layout members */
  layout->latsize         = lattice_size;
  layout->latdim          = LATDIM;
  layout->volume          = volume;
  layout->sites_on_node   = sites_on_node;
  layout->this_node       = this_node;
  layout->number_of_nodes = number_of_nodes;
}

/* Default rank for serial I/O */
static int master_node_zero(){
  return 0;
}

void build_qio_filesystem(QIO_Filesystem *fs){
  fs->number_io_nodes = 0;
  fs->type = QIO_SINGLE_PATH;
  fs->my_io_node = io_node;        /* Partfile I/O uses io_node from layout*.c */
  fs->master_io_node = master_node_zero;  /* Serial I/O uses default: node 0 */
  fs->io_node = NULL;
  fs->node_path = NULL;
}

/* Translate QIO status to a MILC convention */
/* 0 = success
  -1 = end of file
   1 = other failure
*/

int qio_status(int status){
  if(status == QIO_SUCCESS)return 0;
  if(status == QIO_EOF)return -1;
  return 1;
}

QIO_Writer *open_scidac_output(const char *filename, int volfmt, 
			       int serpar, int ildgstyle, 
			       const char *stringLFN, QIO_Layout *layout,
			       QIO_Filesystem *fs,
			       QIO_String *xml_write_file){
  QIO_Writer *outfile;
  QIO_Oflag oflag;

  /* Create the output flag structure */
  oflag.serpar = serpar;
  oflag.ildgstyle = ildgstyle;
  if(stringLFN != NULL){
    oflag.ildgLFN = QIO_string_create();
    QIO_string_set(oflag.ildgLFN, stringLFN);
  }
  else
    oflag.ildgLFN = NULL;
  oflag.mode = QIO_TRUNC;

  /* Open the file for writing */
#ifdef QIO_TRELEASE
  QIO_set_trelease(0,QIO_TRELEASE);
#endif
  QIO_verbose(QIO_VERB_OFF);
  outfile = QIO_open_write(xml_write_file, filename, volfmt, layout, 
			   fs, &oflag);
  if(outfile == NULL){
    printf("open_scidac_output(%d): QIO_open_write returned NULL\n",this_node);
    fflush(stdout);
    return NULL;
  }
  return outfile;
}

/* Open file silently and return file XML */

QIO_Reader *open_scidac_input_xml(const char *filename, QIO_Layout *layout,
				  QIO_Filesystem *fs, int serpar,
				  QIO_String *xml_file_in){
  QIO_Reader *infile;
  QIO_Iflag iflag;
  char myname[] = "open_scidac_input_xml";

  /* Create the iflag structure */
  iflag.serpar = serpar;
  iflag.volfmt = QIO_UNKNOWN;  /* Just discover the format */

  /* Open the file for reading */
#ifdef QIO_TRELEASE
  QIO_set_trelease(0,QIO_TRELEASE);
#endif
  infile = QIO_open_read(xml_file_in, filename, layout, fs, &iflag);
  if(infile == NULL){
    printf("%s(%d): QIO_open_read returns NULL.\n",myname,this_node);
    fflush(stdout);
    return NULL;
  }
  return infile;
}

/* Open file with announcement and discard file XML */

QIO_Reader *open_scidac_input(const char *filename, QIO_Layout *layout, 
			      QIO_Filesystem *fs, int serpar){
  QIO_String *xml_file_in;
  QIO_Reader *infile;

  /* Allocate for the file XML string */
  xml_file_in = QIO_string_create();

  infile = open_scidac_input_xml(filename, layout, fs, serpar, xml_file_in);

  if(infile == NULL)return NULL;

  if(this_node==0){
    printf("Restoring binary SciDAC file %s\n",filename);
    printf("File info \n\"%s\"\n",QIO_string_ptr(xml_file_in));
  }

  QIO_string_destroy(xml_file_in);
  return infile;
}

void close_scidac_output(QIO_Writer *outfile)
{
  QIO_close_write(outfile);
}

void close_scidac_input(QIO_Reader *infile)
{
  QIO_close_read(infile);
}


/********************************************************************/

#ifdef NO_GAUGE_FIELD
static gauge_file *
save_scidac(su3_matrix *field,  const char *filename, int volfmt, int serpar, int ildgstyle,
	    int prec, const char *stringLFN){
  printf("Can't save a lattice if we compile with -DNO_GAUGE_FIELD\n");
  terminate(1);
  return NULL;
}
#else
/* Save the lattice in the site structure in SciDAC format */
/* The QIO file is closed after writing the lattice */
static gauge_file *
save_scidac(su3_matrix *field, const char *filename, int volfmt, int serpar, int ildgstyle,
	    int prec, const char *stringLFN){
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  int status;
  gauge_file *gf;
  char *info;
  QIO_String *filexml;
  QIO_String *recxml;
  char default_file_xml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>MILC ILDG archival gauge configuration</title>";

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O system */
  build_qio_filesystem(&fs);

  /* Make a dummy gauge file structure for MILC use */
  gf = setup_output_gauge_file();

  /* Set the filename in the gauge_file structure */
  gf->filename = filename;

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml, default_file_xml);
  outfile = open_scidac_output(filename, volfmt, serpar, ildgstyle, 
			       stringLFN, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);

  /* Create the QCDML string for this configuration */
  info = create_QCDML();
  recxml = QIO_string_create();
  QIO_string_set(recxml, info);

  /* Write the lattice field */
  if(field == NULL){
    field_offset src = F_OFFSET(link[0]);
    if(prec == 1)
      status = write_F3_M_from_site(outfile, recxml, src, LATDIM);
    else
      status = write_D3_M_from_site(outfile, recxml, src, LATDIM);
  } else {
    if(prec == 1)
      status = write_F3_M_from_field(outfile, recxml, field, LATDIM);
    else
      status = write_D3_M_from_field(outfile, recxml, field, LATDIM);
  }

  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved gauge configuration as a single binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved gauge configuration as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved gauge configuration in partition format to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE_DIR){
    node0_printf("Saved gauge configuration in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Time stamp %s\n",gf->header->time_stamp);
  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  /* Close the file */
  QIO_close_write(outfile);

  free_QCDML(info);
  return gf;
}
#endif /* NO_GAUGE_FIELD */


/* The functions below constitute the API */

int read_lat_dim_scidac(const char *filename, int *ndim, int dims[])
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  int i;
  int *latsize;
  QIO_Reader *infile;

  QIO_verbose(QIO_VERB_REG);

  /* Build the layout structure */
  nx = 0; ny = 0; nz = 0; nt = 0;
  build_qio_layout(&layout);
  /* Forces discovery */
  layout.latdim = 0;

  /* Set the file system parameters */
  build_qio_filesystem(&fs);

  /* Get lattice dimensions from file */
  infile = open_scidac_input(filename, &layout, &fs, QIO_SERIAL);
  if(!infile)return 1;

  *ndim = QIO_get_reader_latdim(infile);
  latsize = QIO_get_reader_latsize(infile);

  for(i = 0; i < *ndim; i++)
    dims[i] = latsize[i];

  QIO_close_read(infile);

  return 0;
}

gauge_file *save_serial_scidac(su3_matrix *field, const char *filename, int prec){
  return save_scidac(field, filename, QIO_SINGLEFILE, QIO_SERIAL, QIO_ILDGNO, prec, NULL);
}

gauge_file *save_parallel_scidac(su3_matrix *field, const char *filename, int prec){
  return save_scidac(field, filename, QIO_SINGLEFILE, QIO_PARALLEL, QIO_ILDGNO, prec, NULL);
}

gauge_file *save_multifile_scidac(su3_matrix *field, const char *filename, int prec){
  return save_scidac(field, filename, QIO_MULTIFILE, QIO_SERIAL, QIO_ILDGNO, prec, NULL);
}

gauge_file *save_partfile_scidac(su3_matrix *field, const char *filename, int prec){
  return save_scidac(field, filename, QIO_PARTFILE, QIO_SERIAL, QIO_ILDGNO, prec, NULL);
}

gauge_file *save_partfile_dir_scidac(su3_matrix *field, const char *filename, int prec){
  return save_scidac(field, filename, QIO_PARTFILE_DIR, QIO_SERIAL, QIO_ILDGNO, prec, NULL);
}

gauge_file *save_serial_ildg(su3_matrix *field, const char *filename, int prec, const char *stringLFN){
  return save_scidac(field, filename, QIO_SINGLEFILE, QIO_SERIAL, QIO_ILDGLAT, prec,
		     stringLFN);
}

gauge_file *save_parallel_ildg(su3_matrix *field, const char *filename, int prec, const char *stringLFN){
  return save_scidac(field, filename, QIO_SINGLEFILE, QIO_PARALLEL, QIO_ILDGLAT, prec,
		     stringLFN);
}

gauge_file *save_multifile_ildg(su3_matrix *field, const char *filename, int prec, const char *stringLFN){
  return save_scidac(field, filename, QIO_MULTIFILE, QIO_SERIAL, QIO_ILDGLAT, prec,
		     stringLFN);
}

gauge_file *save_partfile_ildg(su3_matrix *field, const char *filename, int prec, const char *stringLFN){
  return save_scidac(field, filename, QIO_PARTFILE, QIO_SERIAL, QIO_ILDGLAT, prec,
		     stringLFN);
}

gauge_file *save_partfile_dir_ildg(su3_matrix *field, const char *filename, int prec, const char *stringLFN){
  return save_scidac(field, filename, QIO_PARTFILE_DIR, QIO_SERIAL, QIO_ILDGLAT, prec,
		     stringLFN);
}

/* The QIO file is closed after reading the lattice */
#if 0
static gauge_file *restore_scidac(const const char *filename, su3_matrix *field, int serpar){
  printf("Can't restore a lattice if we compile with -DNO_GAUGE_FIELD\n");
  terminate(1);
  return NULL;
}

#endif

#ifndef NO_GAUGE_FIELD

static gauge_file *restore_scidac(su3_matrix *field, const char *filename, int serpar){
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_RecordInfo recinfo;
  QIO_String *recxml;
  int status;
  int typesize;
  gauge_file *gf;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);
  node0_printf("Calling build_qio_layout with volume %lu\n", layout.volume);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Make a dummy gauge file structure for MILC use */
  gf = setup_input_gauge_file(filename);

  /* Set the filename in the gauge_file structure */
  gf->filename = filename;

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout, &fs, serpar);
  if(infile == NULL)terminate(1);

  /* Check the record type (double or single precision) */
  recxml = QIO_string_create();
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status)terminate(1);
  typesize = QIO_get_typesize(&recinfo);

  /* Read the lattice field as single or double precision according to
     the type size (bytes in a single SU(3) matrix) */
  if(field == NULL){
    field_offset dest = F_OFFSET(link[0]);
    if(typesize == 72)
      status = read_F3_M_to_site(infile, recxml, dest, LATDIM);
    else if (typesize == 144)
      status = read_D3_M_to_site(infile, recxml, dest, LATDIM);
    else
      {
	node0_printf("restore_scidac: Bad typesize %d\n",typesize);
	terminate(1);
      }
  } else {
    if(typesize == 72)
      status = read_F3_M_to_field(infile, recxml, field, LATDIM);
    else if (typesize == 144)
      status = read_D3_M_to_field(infile, recxml, field, LATDIM);
    else
      {
	node0_printf("restore_scidac: Bad typesize %d\n",typesize);
	terminate(1);
      }
  }
  
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);

  /* Close the file */
  QIO_close_read(infile);

  return gf;
}

/* Read gauge field, but don't store the values.  Used
   for verifying checksums */

static gauge_file *file_scan_scidac(const char *filename, int serpar){
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_RecordInfo recinfo;
  QIO_String *recxml;
  int status;
  int typesize;
  su3_matrix *dest = NULL;
  gauge_file *gf;
  int ndim;
  int dims[4];

  /* Read header to get lattice dimensions and close the file */
  read_lat_dim_scidac(filename, &ndim, dims);
  nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
  volume = (size_t) nx*ny*nz*nt;

  /* Finish setting up, now we know the dimensions */

  broadcast_bytes((char *)&nx,sizeof(int));
  broadcast_bytes((char *)&ny,sizeof(int));
  broadcast_bytes((char *)&nz,sizeof(int));
  broadcast_bytes((char *)&nt,sizeof(int));
  
  setup_layout();

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Make a dummy gauge file structure for MILC use */
  gf = setup_input_gauge_file(filename);

  /* Set the filename in the gauge_file structure */
  gf->filename = filename;

  /* Reopen file for reading */
  QIO_verbose(QIO_VERB_OFF);

  infile = open_scidac_input(filename, &layout, &fs, serpar);
  if(infile == NULL)terminate(1);

  /* Check the record type (double or single precision) */
  recxml = QIO_string_create();
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status)terminate(1);
  typesize = QIO_get_typesize(&recinfo);

  /* Read the lattice field as single or double precision according to
     the type size (bytes in a single SU(3) matrix) */
  if(typesize == 72)
    status = read_F3_M_to_null(infile, recxml, dest, LATDIM);
  else if (typesize == 144)
    status = read_D3_M_to_null(infile, recxml, dest, LATDIM);
  else
    {
      node0_printf("file_scan_scidac: Bad typesize %d\n",typesize);
      terminate(1);
    }
  if(status){
    node0_printf("ERROR scanning file\n");
  } else {
    node0_printf("SUCCESS scanning file\n");
  }

  /* Discard for now */
  QIO_string_destroy(recxml);

  /* Close the file */
  QIO_close_read(infile);

  return gf;
}

/* In all cases the QIO file is closed after reading the lattice */
gauge_file *restore_serial_scidac(su3_matrix *field, const char *filename){
  return restore_scidac(field, filename, QIO_SERIAL);
}

gauge_file *restore_parallel_scidac(su3_matrix *field, const char *filename){
  return restore_scidac(field, filename, QIO_PARALLEL);
}

gauge_file *file_scan_serial_scidac(const char *filename){
  return file_scan_scidac(filename, QIO_SERIAL);
}

gauge_file *file_scan_parallel_scidac(const char *filename){
  return file_scan_scidac(filename, QIO_PARALLEL);
}

#endif

/* Read color matrices in SciDAC format */
void restore_color_matrix_scidac_to_site(const char *filename, 
		 field_offset dest, int count){
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *recxml;
  int status, typesize;
  QIO_RecordInfo recinfo;
  char myname[] = "restore_color_matrix_scidac_to_site";

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Set the file system parameters */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout, &fs, QIO_SERIAL);
  if(infile == NULL)terminate(1);

  /* Check the record type (double or single precision) */
  recxml = QIO_string_create();
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status != QIO_SUCCESS)terminate(1);
  typesize = QIO_get_typesize(&recinfo);

  /* Read the lattice field */
  if(typesize == 18*4)
    status = read_F3_M_to_site(infile, recxml, dest, count);
  else if(typesize == 18*8)
    status = read_D3_M_to_site(infile, recxml, dest, count);
  else {
    node0_printf("%s: Incorrect data type size %d\n", myname, typesize);
    terminate(1);
  }
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);

  close_scidac_input(infile);
}

/* Read color matrices in SciDAC format to a field */
void restore_color_matrix_scidac_to_field(const char *filename, 
		  su3_matrix *dest, int count, int prec){
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *recxml;
  int status, typesize;
  QIO_RecordInfo recinfo;
  char myname[] = "restore_color_matrix_scidac_to_field";

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Set the file system parameters */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout, &fs, QIO_SERIAL);
  if(infile == NULL)terminate(1);

  /* Check the record type (double or single precision) */
  recxml = QIO_string_create();
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status != QIO_SUCCESS)terminate(1);
  typesize = QIO_get_typesize(&recinfo);

  /* Read the lattice field */
  if(typesize == 18*4)
    status = read_F3_M_to_field(infile, recxml, dest, count);
  else if(typesize == 18*8)
    status = read_D3_M_to_field(infile, recxml, dest, count);
  else {
    node0_printf("%s: Incorrect data type size %d prec %d\n", myname, typesize, prec);
    terminate(1);
  }
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);

  close_scidac_input(infile);
}

/* Write a set of color matrices in SciDAC format, taking data from the site
   structure */
void save_color_matrix_scidac_from_site(const char *filename, const char *fileinfo, 
		const char *recinfo, int volfmt,  field_offset src, int count, int prec,
		const char *stringLFN)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_String *filexml;
  QIO_String *recxml;
  int ildgType;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O system */
  build_qio_filesystem(&fs);

  /* Is this ILDG format? */
  if(stringLFN != NULL)
    ildgType = QIO_ILDGLAT;
  else
    ildgType = QIO_ILDGNO;

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml,fileinfo);
  outfile = open_scidac_output(filename, volfmt, QIO_SERIAL,
			       ildgType, stringLFN, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);

  /* Write the lattice field */
  recxml = QIO_string_create();
  QIO_string_set(recxml, recinfo);
  if(prec == 1)
    status = write_F3_M_from_site(outfile, recxml, src, count);
  else
    status = write_D3_M_from_site(outfile, recxml, src, count);
  if(status)terminate(1);
  QIO_string_destroy(recxml);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved KS matrix serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved KS matrix as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved KS matrix in partition format to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE_DIR){
    node0_printf("Saved KS matrix in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  close_scidac_output(outfile);
}


/* Save a set of color matrices. */

void save_color_matrix_scidac_from_field(const char *filename,
        const char *fileinfo, const char *recinfo, int volfmt, su3_matrix *src, 
	int count, int prec, const char *stringLFN)
{
  QIO_Layout layout;
  QIO_Writer *outfile;
  QIO_Filesystem fs;
  QIO_String *filexml;
  QIO_String *recxml;
  int ildgType;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Build the structure defining the I/O nodes */
  build_qio_filesystem(&fs);

  /* Is this ILDG format? */
  if(stringLFN != NULL)
    ildgType = QIO_ILDGLAT;
  else
    ildgType = QIO_ILDGNO;

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml, fileinfo);
  outfile = open_scidac_output(filename, volfmt, QIO_SERIAL,
                               ildgType, stringLFN, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);

  /* Write the lattice field */
  recxml = QIO_string_create();
  QIO_string_set(recxml, recinfo);
  if(prec == 1)
    status = write_F3_M_from_field(outfile, recxml, src, count);
  else
    status = write_D3_M_from_field(outfile, recxml, src, count);
  if(status)terminate(1);
  QIO_string_destroy(recxml);

  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved KS matrix serially to binary file %s\n",
                 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved KS matrix as multifile to binary file %s\n",
           filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved KS matrix in partition format to binary file %s\n",
           filename);
  }
  else if(volfmt == QIO_PARTFILE_DIR){
    node0_printf("Saved KS matrix in partition format to binary file %s\n",
           filename);
  }

  node0_printf("Checksums %x %x\n",
               QIO_get_writer_last_checksuma(outfile),
               QIO_get_writer_last_checksumb(outfile));

  close_scidac_output(outfile);
}

/* Read random number state in SciDAC format (SITERAND case only) */
void restore_random_state_scidac_to_site(const char *filename, field_offset dest){
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *recxml;
  int status;

#ifndef SITERAND
  node0_printf("restore_random_state_scidac_to_site: requires SITERAND. Save skipped\n");
  if(1)return;
#endif  

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Set the file system parameters */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout, &fs, QIO_SERIAL);
  if(infile == NULL)terminate(1);

  /* Read the lattice field */
  recxml = QIO_string_create();
  status = read_S_to_site(infile, recxml, dest);
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);

  close_scidac_input(infile);
}

/* Save random number state. (SITERAND case only) */
void save_random_state_scidac_from_site(const char *filename, 
        const char *fileinfo, const char *recinfo, int volfmt, field_offset src)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_String *filexml;
  QIO_String *recxml;
  int status;

#ifndef SITERAND
  node0_printf("save_random_state_scidac_from_site: requires SITERAND. Save skipped\n");
  if(1)return;
#endif  

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O system */
  build_qio_filesystem(&fs);

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml, fileinfo);
  outfile = open_scidac_output(filename, volfmt, QIO_SERIAL,
			       QIO_ILDGNO, NULL, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);

  /* Write the lattice field */
  recxml = QIO_string_create();
  QIO_string_set(recxml, recinfo);
  status = write_S_from_site(outfile, recxml, src);
  QIO_string_destroy(recxml);
  if(status)terminate(1);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved random state serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved random state multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved random state in partition format to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE_DIR){
    node0_printf("Saved random state in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  close_scidac_output(outfile);
}

#if 0
/* Restore a real field */

void restore_real_scidac_to_field(const char *filename, Real *dest, int count){
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *recxml;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Set the file system parameters */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout, &fs, QIO_SERIAL);
  if(infile == NULL)terminate(1);

  /* Read the lattice field */
  recxml = QIO_string_create();
  status = read_F_R_to_field(infile, recxml, dest, count);
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);

  close_scidac_input(infile);
}

#endif

/* Write a complex field */

QIO_Writer *w_open_scidac_file(const char *filename, const char *fileinfo, 
			       int volfmt, int serpar)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_String *filexml;

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml, fileinfo);
  outfile = open_scidac_output(filename, volfmt, serpar, QIO_ILDGNO,
			       NULL, &layout, &fs, filexml);
  QIO_string_destroy(filexml);
  return outfile;
}

int save_complex_scidac(QIO_Writer *outfile, const char *filename, char *recinfo,
			int volfmt, complex *src, int count)
{
  QIO_String *recxml;
  int status;

  recxml = QIO_string_create();
  QIO_string_set(recxml, recinfo);
  status = write_F_C_from_field(outfile, recxml, src, count);
  QIO_string_destroy(recxml);
  if(status)return status;
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved complex field serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved complex field as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved complex field in partition format to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE_DIR){
    node0_printf("Saved complex field in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  return status;
}

void w_close_scidac_file(QIO_Writer *outfile)
{
  QIO_close_write(outfile);
}

/* Save a complex field */

void save_complex_scidac_from_field(const char *filename, char *fileinfo,
				    char *recinfo, int volfmt, int serpar, 
				    complex *src, int count){
    QIO_Writer *outfile;
    int status;

    QIO_verbose(QIO_VERB_OFF);

    outfile = w_open_scidac_file(filename, fileinfo, volfmt, serpar);
    if(outfile == NULL)terminate(1);

    /* Write the lattice field: "count" complex numbers per site */
    status = save_complex_scidac(outfile, filename, recinfo,
				 volfmt, src, count);
    if(status)terminate(1);

    w_close_scidac_file(outfile);
}

/* Restore a complex field */

QIO_Reader *r_open_scidac_file_xml(const char *filename, int serpar,
				   QIO_String *xml_file)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Set the file system parameters */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input_xml(filename, &layout, &fs, serpar, xml_file);

  if(infile == NULL)return NULL;

  if(this_node==0){
    printf("Restoring binary SciDAC file %s\n",filename);
    printf("File info \n\"%s\"\n",QIO_string_ptr(xml_file));
  }

  return infile;
}

QIO_Reader *r_open_scidac_file(const char *filename, int serpar)
{
  QIO_Reader *infile;
  QIO_String *xml_file;

  /* Open file for reading */
  xml_file = QIO_string_create();
  infile = r_open_scidac_file_xml(filename, serpar, xml_file);
  QIO_string_destroy(xml_file);
  return infile;
}

/* Restore a complex field */

int read_complex_scidac_xml(QIO_Reader *infile, complex *dest, int count, 
			    QIO_String *recxml){
  int status;
  int typesize, datacount;
  QIO_RecordInfo recinfo;

  /* Read the lattice field: "count" complex numbers per site */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  status = qio_status(status);
  if(status < 0)return -1;
  if(status > 0)terminate(1);

  datacount = QIO_get_datacount(&recinfo);
  if(datacount != count){
    node0_printf("read_complex_scidac_xml: Got datacount %d but wanted %d\n",
		 datacount, count);
    terminate(1);
  }
  typesize = QIO_get_typesize(&recinfo);

  if(typesize == 2*4)
    status = read_F_C_to_field(infile, recxml, dest, count);
  else
    status = read_D_C_to_field(infile, recxml, dest, count);

  return status;
}

int read_complex_scidac(QIO_Reader *infile, complex *dest, int count){
  QIO_String *recxml;
  int status;

  /* Read the lattice field: "count" complex numbers per site */
  recxml = QIO_string_create();

  status = read_complex_scidac_xml(infile, dest, count, recxml);

  /* Discard for now */
  QIO_string_destroy(recxml);

  return status;
}

void r_close_scidac_file(QIO_Reader *infile)
{
  close_scidac_input(infile);
}

/* Restore a complex field */

void restore_complex_scidac_to_field(const char *filename, int serpar,
				     complex *dest, int count){
  QIO_Reader *infile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  infile = r_open_scidac_file(filename, serpar);
  if(infile == NULL)terminate(1);
  
  /* Read the lattice field: "count" complex numbers per site */
  status = read_complex_scidac(infile, dest, count);
  if(status)terminate(1);

  r_close_scidac_file(infile);
}


/* Restore a real field */

int read_real_scidac_xml(QIO_Reader *infile, Real *dest, int count, 
			    QIO_String *recxml){
  int status;
  int typesize, datacount;
  QIO_RecordInfo recinfo;

  /* Read the lattice field: "count" complex numbers per site */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  status = qio_status(status);
  if(status < 0)return -1;
  if(status > 0)terminate(1);

  datacount = QIO_get_datacount(&recinfo);
  if(datacount != count){
    node0_printf("read_real_scidac_xml: Got datacount %d but wanted %d\n",
		 datacount, count);
    terminate(1);
  }
  typesize = QIO_get_typesize(&recinfo);

  if(typesize == 4)
    status = read_F_R_to_field(infile, recxml, dest, count);
  else
    status = read_D_R_to_field(infile, recxml, dest, count);

  return status;
}

int read_real_scidac(QIO_Reader *infile, Real *dest, int count){
  QIO_String *recxml;
  int status;

  /* Read the lattice field: "count" complex numbers per site */
  recxml = QIO_string_create();

  status = read_real_scidac_xml(infile, dest, count, recxml);

  /* Discard for now */
  QIO_string_destroy(recxml);

  return status;
}

/* Restore a real field */

void restore_real_scidac_to_field(const char *filename, int serpar,
				  Real *dest, int count){
  QIO_Reader *infile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  infile = r_open_scidac_file(filename, serpar);
  if(infile == NULL)terminate(1);
  
  /* Read the lattice field: "count" complex numbers per site */
  status = read_real_scidac(infile, dest, count);
  if(status)terminate(1);

  r_close_scidac_file(infile);
}


/* Restore a real field */

void restore_real_scidac_to_site(const char *filename, field_offset dest, int count)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *recxml;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Set the file system parameters */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout, &fs, QIO_SERIAL);
  if(infile == NULL)terminate(1);

  /* Read the lattice field */
  recxml = QIO_string_create();
  status = read_F_R_to_site(infile, recxml, dest, count);
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);

  close_scidac_input(infile);
}


/* Save a real field */

void save_real_scidac_from_field(const char *filename, 
        const char *fileinfo, const char *recinfo, int volfmt, Real *src, int count)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_String *filexml;
  QIO_String *recxml;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O system */
  build_qio_filesystem(&fs);

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml,fileinfo);
  outfile = open_scidac_output(filename, volfmt, QIO_SERIAL,
			       QIO_ILDGNO, NULL, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);

  /* Write the lattice field */
  recxml = QIO_string_create();
  QIO_string_set(recxml,recinfo);
  status = write_F_R_from_field(outfile, recxml, src, count);
  if(status)terminate(1);
  QIO_string_destroy(recxml);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved real field serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved real field in multifile format to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved real field in partition format to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE_DIR){
    node0_printf("Saved real field in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  close_scidac_output(outfile);
}


/* Save a real field from the MILC site structure*/

void save_real_scidac_from_site(const char *filename, 
  const char *fileinfo, const char *recinfo, int volfmt, field_offset src, int count)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_String *filexml;
  QIO_String *recxml;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O system */
  build_qio_filesystem(&fs);

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml,fileinfo);
  outfile = open_scidac_output(filename, volfmt, QIO_SERIAL,
			       QIO_ILDGNO, NULL, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);

  /* Write the lattice field */
  recxml = QIO_string_create();
  QIO_string_set(recxml,recinfo);
  status = write_F_R_from_site(outfile, recxml, src, count);
  if(status)terminate(1);
  QIO_string_destroy(recxml);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved real field serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved real field in multifile format to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved real field in partition format to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE_DIR){
    node0_printf("Saved real field in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  close_scidac_output(outfile);
}



