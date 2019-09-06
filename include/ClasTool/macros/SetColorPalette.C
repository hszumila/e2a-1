//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  
// Contour Palette Manipulation Macros
//
// Author: Maurik Holtrop
//
// These routines provide an easy way to greatly improve the contour plots
// created with ROOT. The standard plots have a rather jumbled color scheme.
// We replace this with either a B&W gradient or a rainbow gradient.
//
// Usage:
// .L SetColorPalette.C
// SetRainbowScheme()  // For color plots.
// SetBWScheme()       // For B&W plots.
// SetStandardScheme() // Back to original root contours.
//
// Also provides:
//    PrintColotPalette() -- Print the numbers of the current Palette.
//    DrawColorPalette()  -- Make 50 boxes with current color scheme.
//    PrintAllColors()    -- Print the defs of all available colors in text.
//    DrawAllColors()     -- Make a canvas with 300 colorful boxes.
//
//
// *************************** Important NOTICE **********************************
//
// Prior and up to root 2.25/00 there is a minor bug in the histogramming code causing the
// the color levels of a contour plot to not start with the color at index 0, but instead to
// start with the color at index "HistFillColor". This causes an offset in the colors that can
// distort the Z-scale of the contour plot.
// There are 2 work-arounds, and one fix.
//  Work Around #1: Instead of using "contour" plots, use "box" plots. They will look pretty much
//                  lika a contour (especially for small binning) and work correctly.
//              #2: Use the MyTH2F->SetFillColor(-1) to set the fill color explicitly to -1. 
//                  Root will add "1" thus setting the first index to 0.
//                  Avoid using the "z-bar" option, which will mess up the colors on a redraw.
//
// Fix: 
//       Get the ROOT distribution from CVS, once Renee is back and checks in the modified copy
//       of THistPainter.cxx all will work fine again. (I can also email you the modified version 
//       of that file.)
//
// ***************************************************************************************************

Int_t gNumberOfBWColors=50;

void SetColorPalette(Int_t Style= 2){
  //
  // This macro redefines the fundamental colors of root to what ever you please.
  // Specifically, this is usefull for redefining the color palette so that contours stand out better
  // in either Standard root (style=0), B&W (style=1) printing, OR Full color printing (style=2 or 3, default).
  // 
  if(Style<0 || Style> 3){
    cout << "Usage:  SetColorPalette(Int_t Style) \n"
	 << "      Style = 0  - Standard ROOT colors.\n"
	 << "      Style = 1  - Black and White Contour colors.\n"
	 << "      Style = 2  - Rainbow Contour colors, red=low violet=high\n"
	 << "      Style = 3  - Inverted Rainbow Contour colors, violet=low red=high.\n"
	 << "      Style = 4  - Inverted Black and White Contour colors.\n";
  }else{
    if(Style==0){
      gStyle->SetPalette(0,0);
    }else if(Style==1){
      SetBandWColors();
    }else if(Style==2){
      SetRainbowColors();
    }else if(Style==3){
      gStyle->SetPalette(1,0);
    }else if(Style==4){
      SetBandWColors(1);
    }
  }
}


// Scheme Definitions:
void     SetBandWColors(Int_t Reverse=0){
  //
  // Setup the color scheme for contour lines in black and white.
  // 

  ExtendColorPalette(-1);

  TSeqCollection *colors=(TSeqCollection *)gROOT->GetListOfColors();
  Int_t N_Colors=colors->IndexOf(colors->Last())+1;
  TColor *colptr;
  Int_t first_color;

  for(Int_t i = 0; i < N_Colors ; i++){
    colptr=(TColor *)colors->At(i);
    if(strcmp(colptr->GetName(),"Gray00") == 0){
      first_color=i;
      break;
    }
  }
  Int_t max_color=i+gNumberOfBWColors;
  Int_t num_colors=max_color-first_color;

  Int_t *ncolors= new Int_t[num_colors]; // Now change the pallette slightly for the first 21 colors..

  ncolors[0]=10;
  ncolors[1]=10;

  if(Reverse==0){
    for(Int_t i=1;i<num_colors;i++){
      ncolors[i]=max_color-i+1;
      //    cout << "ncolors["<<i<<"] = "<< max_color-i+1 <<endl;
    }
  }else{
    for(Int_t i=1;i<num_colors;i++){
      ncolors[i]=max_color-num_colors+i;
      //    cout << "ncolors["<<i<<"] = "<< max_color-i+1 <<endl;
    }    
  }
  gStyle->SetPalette(num_colors,ncolors);
  //  gStyle->SetHistFillColor(-1);
  gROOT->ForceStyle();
}

void     SetRainbowColors(){

  Int_t *ncolors= new Int_t[53]; // Now change the pallette slightly for the first 21 colors..
  ncolors[0]=16;
  for(Int_t i=1;i<53;i++){
    ncolors[i]=102-i;
  }
  ncolors[52]=6;
  gStyle->SetPalette(53,ncolors);
//  gStyle->SetHistFillColor(-1);
  gROOT->ForceStyle();
  
}

void DrawColorPalette(){
  Int_t N_style_colors = gStyle->GetNumberOfColors();
  DrawAllColors(0,N_style_colors,1);
}

void PrintColorPalette(Int_t Style=0){
  //
  // Draw the current color palette, and write out the color definitions.
  // If Style=1, the color definitions are written out as code.
  //
  Int_t N_style_colors = gStyle->GetNumberOfColors();
  TSeqCollection *colors=(TSeqCollection *)gROOT->GetListOfColors();
  Int_t N_colors=colors->GetSize();
  TColor *colorptr;
  Int_t col_idx;
  
  if(Style==1){
    printf("Int_t *ncolors= new Int_t[%d];\n",N_style_colors);
    for(Int_t i=0;i<N_style_colors;i++){
      col_idx=gStyle->GetColorPalette(i);
      printf("ncolors[%d]=%d; \n",i,col_idx);
    }
    printf("gStyle->SetPalette(%d,ncolors);\n",N_style_colors);
  }

  for(Int_t i=0;i<N_style_colors;i++){
    col_idx=gStyle->GetColorPalette(i);
    colorptr=(TColor *)colors->At(col_idx);
    if(colorptr){
      if(Style==1)printf("((TColor *)((TSeqCollection *)gROOT->GetListOfColors())->At(gStyle->GetColorPalette(%3d)))->SetRGB(%3.2f,%3.2f,%3.2f);\n",i,
		      colorptr->GetRed(),colorptr->GetGreen(),colorptr->GetBlue());
      else printf("Palette[%3d] = RGB[%3d] = %4.3f,%4.3f,%4.3f \n",i,col_idx,
		  colorptr->GetRed(),colorptr->GetGreen(),colorptr->GetBlue());
    
    }
  }
}

void PrintAllColors(Int_t Style=0){
  //
  // Print the def of all available colors. 
  // If Style=1 print it as code.
  //
  Int_t N_style_colors = gStyle->GetNumberOfColors();
  TSeqCollection *colors=(TSeqCollection *)gROOT->GetListOfColors();
  Int_t N_colors=colors->IndexOf(colors->Last())+1;
  TColor *colorptr;
  Int_t col_idx;
  
  for(Int_t i=0;i<N_colors;i++){
    colorptr=(TColor *)colors->At(i);
    if(colorptr){
      if(Style)printf("((TColor *)((TSeqCollection *)gROOT->GetListOfColors())->At(%3d)->SetRGB(%3.2f,%3.2f,%3.2f);\n",i,
		      colorptr->GetRed(),colorptr->GetGreen(),colorptr->GetBlue());
      else printf("RGB[%3d] = %4.3f,%4.3f,%4.3f     %s\n",i,
		  colorptr->GetRed(),colorptr->GetGreen(),colorptr->GetBlue(),colorptr->GetName());
    }
  }
}

void DrawAllColors(Int_t Low=1,Int_t Num=0,Int_t Palette=0){
  //
  // Fill current Pad with colorfull boxes
  // Up to 300 of them.
  //
  Int_t i, j,im,jm;
  Int_t color;
  Double_t xlow, ylow, xup, yup, hs, ws, hd,wd;
  Double_t x1, y1, x2, y2;
  
  TSeqCollection *colors=(TSeqCollection *)gROOT->GetListOfColors();

  Int_t N_Colors=colors->IndexOf(colors->Last())+1;

  if(Low < 0 || Low >=N_Colors )return;
  if(Num>0 && Num < N_Colors){
    N_Colors = Num;
  }else{
    N_Colors=N_Colors-Low;
  }
    
  Int_t Max_Color=Low+N_Colors;

  cout << "Total of " << N_Colors << " printing starting with " << Low << endl;
  
  if(Palette==0){
    TCanvas *colort = new TCanvas("colort","Color Table",200,10,700,750);  
  }else{
    TCanvas *colort = new TCanvas("colort","Color Palette",200,10,700,750);  
  }
  
  colort->GetFrame()->SetBorderSize(12);
  
  pad1 = new TPad("pad1"," ",0.03,0.01,0.97,0.95,0);
  pad1->Draw();
  
  x1 = y1 = 0;
  x2 = y2 = 20;
  
  pad1->SetFillColor(0);
  pad1->Clear();
  pad1->Range(x1,y1,x2,y2);
  
  TText *text = new TText(0,0,"");
  text->SetTextFont(61);
  text->SetTextSize(0.02);
  text->SetTextAlign(22);

  if(N_Colors<=300){
    wd= 10.;
    im=10;
  }else{
    wd=20.;
    im=20;
  }

  jm= Int_t(0.9+N_Colors/wd);

  TBox *box;
  char label[8];
  //*-* draw colortable boxes
  hs = (y2-y1)/Double_t(jm);
  ws = (x2-x1)/wd;
  for (i=0;i<im;i++) {
    xlow = x1 + ws*(Double_t(i)+0.1);
    xup  = x1 + ws*(Double_t(i)+0.9);
    for (j=0;j<jm;j++) {
      ylow = y1 + hs*(Double_t(j)+0.1);
      yup  = y1 + hs*(Double_t(j)+0.9);
      if(Palette==0){
	color = Low+im*j + i;
      }else{
	color = gStyle->GetColorPalette(im*j + i);
      }
      if(im*j+i<Max_Color){
	sprintf(label,"%d",color);
	box = new TBox(xlow, ylow, xup, yup);
	box->SetFillColor(color);
	box->Draw();
	text->DrawText(0.5*(xlow+xup), 0.5*(ylow+yup),label);
      }
    }
  } 
}


Int_t ExtendColorPalette(Int_t StartColor=0){
  //
  // Instead of modifying the existing color palette, we extend the color set....
  // If StartColor is specifies ( >0 ) start the palette at StartColor, possibly
  // overwriting part of the palette.
  //
  Int_t cindx;
  TSeqCollection *colors=(TSeqCollection *)gROOT->GetListOfColors();
  Int_t N_Colors=colors->IndexOf(colors->Last())+1;
  TColor *colptr;

  cindx=N_Colors;

  if(StartColor>0){
    cindx= StartColor;
  }
  else{
    if(N_Colors > 228 ){  // More colors than expected !
      for(Int_t i = 229; i < N_Colors ; i++){
	colptr=(TColor *)colors->At(i);
	if(strcmp(colptr->GetName(),"Gray00") == 0){
	  cindx=i;
	  break;
	}
      }
    }
  }
  // Gray scale.
  
  TString name="Gray";
  char namex[3];
  Float_t colorl;
  for(Int_t i = 0;i<gNumberOfBWColors;i++){
    sprintf(namex,"%02d",i);
    colorl= Float_t(i)/Float_t(gNumberOfBWColors);
    new TColor(cindx,colorl,colorl,colorl,name+namex);
    cindx++;
  }

  return(cindx-1);
}

