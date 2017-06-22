# the original CampSpec IRAF script by C Kulesa et al

procedure campspec (images, output)

string    images      {prompt="name of input image"}
string    output      {prompt="name of output spectrum"}
bool      darksub     {prompt="Perform dark-current subtraction?"}
string    darkimg     {prompt="Filename of dark-current image"}
bool      displayvar  {prompt="Display image frames while processing?"}
string    linelist    {"/home/ckulesa/scripts/HgNe.dat",prompt="Where's your Hg/Ne linelist?"}
struct    *flist

begin

   string outfile, image,drk,list
   int i,n,m
   bool msky,shw,dsp,ds
   real lsig,hsig,lthresh,hthresh

   set imtype = "fits"

# Copy task parameters to local variable names so they will only be queried
# for once.
   outfile=output
   image=images
   ds=darksub
   drk=darkimg
   dsp=displayvar
   list=linelist

# Check to see if images exist
   if(!access(image))
      error(1,"Help! I cannot find ("//image//").")
   if(ds)
   {
      if(!access(drk))
          error(1,"Help! I cannot find your dark current image!")
   }
# Check to see if output file already exists
   if(access(outfile)||access(outfile//'.imh'))
      error(1,"Uh oh! The file ("//outfile//") already exists!")

# The following three lines allow wildcards, @lists and comma-delimited image
# lists to be used - consistent with regular IRAF tasks
#   namelist=mktemp('tmp$tmp')
#   files(images,sort+,>namelist)
#   flist=namelist
   
   if(access("tmplamps.fits"))
   {
      print("Cleaning up temporary files from a previous run...")
      imdel("tmplamps")
   }
   if(access("tmpextract.fits"))
   {
      print("Cleaning up temporary files from a previous run...")
      imdel("tmpextract")
   }
   if(dsp)
   {
      print("Displaying original image...")
      display(image=image,frame=1,zscale=yes,zrange=yes,contras=0.25)
   }
   if(ds)
   {
      print("Now subtracting a dark-current image from your spectrum")
      imarith(operand1=image,op='-',operand2=drk,result=image,verbose=yes)
      if(dsp)
      {
          print("Dark subtraction done. I'm displaying the result now.")
          display(image=image,frame=1,zscale=yes,zrange=yes,contras=0.25)
      }
   }   
   print("Now it's time to extract lamp and object spectra from your two-")
   print("dimensional image. Extract the lamps in aperture 1, and the object")
   print("in aperture 2.")
      
   apall(input=image,output = 'tmpextract',inter=yes,find=no,
      recenter=no,resize=no,edit=yes,trace=yes ,fittrac=yes,
      extract=yes,extras=no,review=no)
   imcopy(input = 'tmpextract[*,1]',output = 'tmplamps')
   identify(images = 'tmplamps',coordli = list)
   reid(referenc = 'tmplamps',images = 'tmplamps',interac=no,refit=no)
   refspec(input = 'tmpextract',referen='tmplamps',sort='',group='',
      confirm=no)
   dispcor(input = 'tmpextract', output=outfile,
      table = list,ignorea=yes,flux=yes,
      samedis=no)
   splot(images=outfile)
   imdel("tmpextract")
   imdel("tmplamps")
   flist = ""

end
