package base_RayTracer.utils;

import java.io.File;
import java.util.*;

import base_Render_Interface.IRenderInterface;
import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.doubles.myVector;
import base_RayTracer.scene.*;
import base_RayTracer.scene.base.Base_Scene;
import base_RayTracer.scene.geometry.sceneObjects.myRndrdBox;
import base_RayTracer.scene.geometry.sceneObjects.base.Base_SceneObject;
import base_RayTracer.scene.geometry.sceneObjects.implicit.*;
import base_RayTracer.scene.geometry.sceneObjects.planar.*;
import base_RayTracer.scene.geometry.sceneObjects.planar.base.Base_PlanarObject;
import base_RayTracer.scene.materials.textures.imageTextures.myImageTexture;
import base_RayTracer.ui.base.Base_RayTracerWin;
import base_UI_Objects.my_procApplet;
import base_UI_Objects.windowUI.base.Base_DispWindow;


public class myRTFileReader {
	public IRenderInterface pa;
	private Base_RayTracerWin win;
	public final String textureDir;
	private int timeSinceStart;
	private int curNumRows, curNumCols;
		
	public myRTFileReader(Base_RayTracerWin _win, String _txtrDirs) {
		win=_win;
		pa = Base_RayTracerWin.ri;
		textureDir = _txtrDirs;		
	}
	
	/**
	 * passed scene is for when called recursively - is null on first time in, passes scene to be used otherwise
	 * @param loadedScenes
	 * @param filePath
	 * @param fileName
	 * @param _scene
	 * @param _numCols
	 * @param _numRows
	 * @return
	 */
	public Base_Scene readRTFile(TreeMap<String, Base_Scene> loadedScenes, String filePath, String fileName, Base_Scene _scene, int _numCols, int _numRows) {		  
		//build individual scene for each file		
		timeSinceStart = getTime();			//times rendering
		curNumRows = _numRows;
		curNumCols = _numCols;
		
		Base_Scene scene = _scene;	
		boolean isMainFileScene = _scene == null;			//whether this is primary scene's file or secondary/recursive read call; whether scene has been named and finalized yet
		if(isMainFileScene){
			//first see if scene exists and has been loaded already (don't reload/remake)
			scene = loadedScenes.get(fileName);
			if(scene != null){
				//if scene exists but is not the same size, do not reload, just re-render with new size
				if(!scene.isSameSize(_numCols, _numRows)) {scene.setNewSize(curNumCols, curNumRows);}				
				return scene;
			}
			//by here we still have null scene (scene was not found in loaded scene map), we load it after we load the file
		}
		String[] strAra = null;
		try{
			//Load the file
			strAra = win.loadRTStrings(filePath + File.separator+fileName);
			win.getMsgObj().dispInfoMessage("myRTFileReader", "readRTFile", "File path : "+filePath + "|File name : " + fileName + " Size : " + strAra.length);
		} catch (Exception e) {	win.getMsgObj().dispErrorMessage("myRTFileReader", "readRTFile", "File Read Error : File name : " + fileName + " not found."); return null;		}	 
		
		//if we make it this far and this is the main scene, read through strAra to build scene.		
		if(isMainFileScene){
			//build the scene described in the cli file
			scene = buildBaseScene(strAra, fileName);
		}
		//populate the scene with the appropriate data from the scene list		
		scene = parseStringArray(loadedScenes, strAra, scene, isMainFileScene, filePath, fileName);
		return scene;
	}//interpreter method
	
	/**
	 * Build the scene described in the .cli file.  
	 * @param fileStrings
	 * @param fileName
	 * @return
	 */
	private Base_Scene buildBaseScene(String[] fileStrings, String fileName) {
		//default scene, if not specified in file
		Base_Scene scene = new myFOVScene(pa, win, fileName, curNumCols, curNumRows, 60.0);
		if (Base_RayTracerWin.AppMgr.isDebugMode()) {	_debugFileStrings(fileStrings);	}
		boolean done = false;
		String sceneType = "fov";
		for(int i=0;i<fileStrings.length;++i){ 
			//Skip comments and empty lines
			if((fileStrings[i].startsWith("#")) || (fileStrings[i].strip().length() == 0)) {continue;}
			String[] tokenAra = fileStrings[i].strip().split("\\s+"); // Get a line and parse tokens.
			//debug
			//_printTokenAra(fileStrings[i], i, tokenAra);
			if ((tokenAra.length == 0) || (tokenAra[0] == "")) {continue;} // Skip blank line or comments.
			switch (tokenAra[0].toLowerCase()){
				//determine the kind of scene
				case "fov" 	: {
					sceneType = tokenAra[0];
					scene = new myFOVScene(pa, win, fileName, curNumCols, curNumRows, Double.parseDouble(tokenAra[1]));
					done = true;
					break;
				}
				case "fisheye" : {
					sceneType = tokenAra[0];
					scene = new myFishEyeScene(pa, win, fileName, curNumCols, curNumRows, Double.parseDouble(tokenAra[1]));
					done = true;
					break;
				}				
				case "ortho" :
				case "orthographic" : {
					sceneType = tokenAra[0];
					scene = new myOrthoScene(pa, win, fileName, curNumCols, curNumRows, Double.parseDouble(tokenAra[1]),Double.parseDouble(tokenAra[2]));
					done = true;
					break;
				}
			}//build scene
			if (done) {
				break;
			}
		}//while not end of list and not done
		if(!done) {
			win.getMsgObj().dispWarningMessage("myRTFileReader", "buildBaseScene", "Warning! No scene type specified in scene file `" + fileName + "` so defaulting to FOV scene with FOV == "+((myFOVScene) scene).getFOV());
		} else {
			win.getMsgObj().dispInfoMessage("myRTFileReader", "buildBaseScene", "Built " + sceneType +" scene for file `" + fileName + "`");
			
		}
		return scene;
	}
	
	/**
	 * TODO use window time manager
	 * @return
	 */
	protected int getTime() {return Base_DispWindow.AppMgr.timeSinceStart();}
	
	/**
	 * build the objects in a scene
	 * @param loadedScenes
	 * @param fileStrings
	 * @param scene
	 * @param isMainFileScene
	 * @param fileName
	 * @return
	 */
	public Base_Scene parseStringArray(TreeMap<String, Base_Scene> loadedScenes, String[] fileStrings, Base_Scene scene, boolean isMainFileScene, String filePath, String fileName) {
		boolean finalized = false;
		String vertType = "triangle";    		//assume default object type is triangle
		int myVertCount = 0;		
		//temp objects intended to hold 
		Base_SceneObject myPoly = null;
		//int curNumRaysPerPxl = scene.numRaysPerPixel;
		//reinitializes the image so that any previous values from other images are not saved
		//if (str == null) {win.getMsgObj().dispErrorMessage("myRTFileReader", "parseStringArray", "Error! Failed to read the file.");}
		//debug display file
		if (Base_RayTracerWin.AppMgr.isDebugMode()) {	_debugFileStrings(fileStrings);	}
		for (int i=0; i<fileStrings.length; ++i) { 
			//Skip comments and empty lines
			if((fileStrings[i].startsWith("#")) || (fileStrings[i].strip().length() == 0)) {continue;}
			String[] tokenAra = fileStrings[i].strip().split("\\s+"); // Get a line and parse tokens.
			//debug
			//_printTokenAra(fileStrings[i], i, tokenAra);
			if ((tokenAra.length == 0) || (tokenAra[0] == "")) {continue;} // Skip blank line or comments.
			switch (tokenAra[0].toLowerCase()){
				//skip scene specifications (already processed)
				case "fov":
				case "fisheye":
				case "ortho":
				case "orthographic": {break;}	//scene already built before this method was called, ignore these elements in file string array
				
				//for depth of field -only in FOV scenes - specifies the lens size (radius) and what 
				//distance in front of the eye is in focus - greater radius should blur more of the image
			    case "lens" : { 			
			    	if(!isMainFileScene){win.getMsgObj().dispErrorMessage("myRTFileReader", "parseStringArray", "Error - unsupported setting scene value ('"+tokenAra[0]+"') in recursive child scene file"); break;}
					double radius = Double.parseDouble(tokenAra[1]);
			    	double focal_distance = Double.parseDouble(tokenAra[2]);
			    	((myFOVScene)scene).setDpthOfFld(radius, focal_distance);			    	
			    	break;}				
				
				//file IO : save and read subordinate file
				//NEEDS TO RENDER SCENE and then save it so that rendering can be timed - doesn't work with refine currently
			    case "write" : {		
			    	scene.setSaveName(tokenAra[1]);			    	
			    	finalizeScene(loadedScenes,fileName, scene);		
			    	finalized = true;
			    	//render scene -here- - needs to be modified to not stop until scene is finished rendering (for incremental/refine scene renders)
			    	if(!fileName.equals("")){scene.draw();}
			    	else{win.getMsgObj().dispErrorMessage("myRTFileReader", "parseStringArray", "Can't render unknown/incomplete scene with empty fileName");}
			    	break;}		
			    case "read" : {//read another scene file - nested to load multiple files
			    	readRTFile(loadedScenes, filePath, tokenAra[1], scene, curNumCols, curNumRows);
			    	break;}					    
			    
			    //timer stuff
			    case "reset_timer" : {//Reset a global timer that will be used to determine how long it takes to render a scene. 
			    	timeSinceStart = getTime();
			    	break;
			    }
			    case "print_timer" : {//Print the time elapsed since reset_timer was called (in seconds). Use the code snippet 
			    	int new_timer = getTime();
			    	int diff = new_timer - timeSinceStart;
			    	float seconds = diff / 1000.0f;
			    	scene.renderTime = seconds;
			    	win.getMsgObj().dispInfoMessage("myRTFileReader", "parseStringArray", "Timer = " + seconds);
			    	break;
			    }			    
			    
			    //global modifications to alg
			    case "refine" : {    	scene.setRefine(tokenAra[1]); break;}//user iterative refinement when rendering scene - start at (int)log2 of dim, then decrease by 2 every iteration			    

			    case "rays_per_pixel" : {	//how many rays should be shot per pixel.
			    	int rays = Integer.parseInt(tokenAra[1]);
			    	win.getMsgObj().dispInfoMessage("myRTFileReader", "parseStringArray","Num Rays Per Pixel : " + rays);
			    	scene.setNumRaysPerPxl(rays);			    	
			    	break;}			
			    
				case "antialias" : {
					int aaDepthRow = Integer.parseInt(tokenAra[1]);
					int aaDepthCol = Integer.parseInt(tokenAra[2]);
					int prod = aaDepthRow * aaDepthCol;
					win.getMsgObj().dispInfoMessage("myRTFileReader", "parseStringArray", "Anti Alias depth r/c " + aaDepthRow + "|" + aaDepthCol +" -> convert to numRays per pixel " + prod);
			    	scene.setNumRaysPerPxl(prod);			    	
					break;} 	
				
				//background texture/color/foreground color
				case "background" : {//visible background in front of camera - negative infinite z\
					if (tokenAra[1].equals("texture")) {
						//load texture to be used for background
						String textureName = tokenAra[2];
						scene.currBkgTexture = ((my_procApplet)pa).loadImage(textureDir+textureName);
						scene.setHasGlblTxtrdBkg(true);
						win.getMsgObj().dispInfoMessage("myRTFileReader", "parseStringArray", "Background texture loaded");
						//build "skydome" - textured sphere encircling scene
						double rad = Double.parseDouble(tokenAra[3]);
						double xC = Double.parseDouble(tokenAra[4]);
						double yC = Double.parseDouble(tokenAra[5]);
						double zC = Double.parseDouble(tokenAra[6]);
						win.getMsgObj().dispInfoMessage("myRTFileReader", "parseStringArray", "Skydome : rad:" + rad +" ctr : ["+ xC +","+ yC +","+ zC+"]");
						scene.mySkyDome = new mySphere(scene,rad,xC,yC,zC,true);
						((myImageTexture)scene.mySkyDome.shdr.txtr).setMyTextureBottom(scene.currBkgTexture);	
					} else {//set bkg color
						scene.setBackgroundColor(Double.parseDouble(tokenAra[1]),Double.parseDouble(tokenAra[2]),Double.parseDouble(tokenAra[3]));
						scene.txtrType = 0;
					}
					break;}
			
				//lights
				case "point_light" : {
					scene.addMyPointLight(tokenAra);	      
					break;}
				case "spotlight" : {//spotlight x y z dx dy dz angle_inner angle_outer r g b
					win.getMsgObj().dispInfoMessage("myRTFileReader", "parseStringArray", "File : " + fileName+ " has a spotLight!");
					scene.addMySpotLight(tokenAra);   
					break;}				
				case "disk_light" : {//disk_light x y z radius dx dy dz r g b
					scene.addMyDiskLight(tokenAra);
					break;}	
				
				//photon mapping
				case "caustic_photons" : //caustic_photons num_cast num_near max_near_dist
				case "diffuse_photons" : {//diffuse_photons num_cast num_near max_near_dist 
					scene.setPhotonHandling(tokenAra);
					break;}
				//
				//final_gather num_rays
				//This command indicates how the diffuse photons for global illumination will be used to create 
				//the final image. If num_rays is set to zero, then your renderer should directly use the diffuse 
				//photons stored on each surface, similar to how you render using caustic photons. 
				//This should be fairly fast, but unfortunately this will create noisy images. If num_rays is non-zero, 
				//then you will estimate the indirect illumination using a "final gather" step. To calculate the indirect 
				//illumination at a surface point, you will shoot num_rays rays in random directions on the hemisphere 
				//surrounding this point. You will then query the kD-tree at each surface that such a ray hit, and use 
				//this to determine how much light should reach the point in question from this surface. This will produce 
				//much better images, but will also be significantly slower.
				case "final_gather" : {//final_gather num_rays
					scene.setFinalGather(tokenAra);
					break;}
			
				//material/color commands
				case "diffuse" : {//new assignment requirements
					myRTColor cDiff = readColor(tokenAra,1);//new myColor(Double.parseDouble(token[1]),Double.parseDouble(token[2]),Double.parseDouble(token[3]));
					myRTColor cAmb = readColor(tokenAra,4);//new myColor(Double.parseDouble(token[4]),Double.parseDouble(token[5]),Double.parseDouble(token[6]));
					myRTColor cSpec = new myRTColor(0,0,0);
					scene.setHasGlblTxtrdTop(false);
					scene.setHasGlblTxtrdBtm(false);
					scene.setSurface(cDiff,cAmb,cSpec,0,0,0);					
					break;}
				//use shiny for new refr/refl; use surface for older cli files - handles mix of colors and refr/refl better
				case "shiny" :{//Cdr Cdg Cdb, Car Cag Cab, Csr Csg Csb, Exp Krefl Ktrans Index 
					setSurfaceShiny(scene, tokenAra, true);
					break;}
				case "surface" : {
					setSurfaceShiny(scene, tokenAra, false);
					break;}
				//reflective Cdr Cdg Cdb Car Cag Cab k_refl 
				//This command describes a new kind of surface material. Just as with the 
				//diffuse command, this command defines the diffuse and ambient color components of a surface. 
				//In addition, however, it also defines a reflectance coefficient k_refl. This value should be in the 
				//range of zero to one. When it is non-zero, k_refl indicates how much of the light that strikes the 
				//surface is to be reflected. Eye rays that hit such a reflective surface should spawn a new reflective ray, 
				//and this reflective ray should contribute to the surface's color. Moreover, if a caustic photon hits such 
				//a reflective surface, this should cause the caustic photon to "bounce", and to travel in a new direction. 
				//Caustic photons only stop bouncing when they hit a diffuse surface. 				
				case "reflective" : {//reflective Cdr Cdg Cdb Car Cag Cab k_refl 
					myRTColor cDiff = readColor(tokenAra,1);
					myRTColor cAmb = readColor(tokenAra,4);
					myRTColor cSpec = new myRTColor(0,0,0);
					scene.setHasGlblTxtrdTop(false);
					scene.setHasGlblTxtrdBtm(false);
					double kRefl = Double.parseDouble(tokenAra[7]);
					scene.setSurface(cDiff,cAmb,cSpec,0,kRefl,0);		
					break;}				
			    case "phong" : {scene.setPhong(Double.parseDouble(tokenAra[1]));	break;}//phong
			    case "perm" : {//load new permiability value - used for refraction			    	
			    	double rfrIdx = Double.parseDouble(tokenAra[1]);
			    	myRTColor permClr = new myRTColor(rfrIdx,rfrIdx,rfrIdx);
			    	if (tokenAra.length>2) {   		permClr = readColor(tokenAra,2);} 
		    		scene.setRfrIdx(rfrIdx,permClr);
			    	break;}//if perm
			    case "krefl" : {
			    	double kRefl = Double.parseDouble(tokenAra[1]);
			    	myRTColor reflClr = new myRTColor(kRefl,kRefl,kRefl);
			    	if (tokenAra.length>2) {   		reflClr = readColor(tokenAra,2);}		    	
			    	scene.setKRefl(kRefl,reflClr);  	
			    	break;}//krefl
				case "depth" : {scene.setDepth(Double.parseDouble(tokenAra[1]));	break;}//depth for subsurface scattering - TODO
			    case "ktrans" : {scene.setKTrans(Double.parseDouble(tokenAra[1]));	break;}//ktrans - //load new permiability value - used for refraction		
			    
			    //accel structs
			    case "begin_list" : {	    	scene.startTmpObjList();		    break;}
			    case "end_list" : {		    	scene.endTmpObjList(0);		    	break;}
			    case "end_accel" : {		   	scene.endTmpObjList(1);		    	break;}			//TODO modify to accept multiple accel struct types		 
			    		    
			    //predefined object layouts
			    case "sierpinski" :{
			    	//generate sierpinski tet arrangement of named object
			    	String objName = tokenAra[1], useShdr="No";
			    	float scale = .5f;
			    	int depth = 5;		//default to 5 layers deep	    	
			    	try {depth = Integer.parseInt(tokenAra[2]);scale = Float.parseFloat(tokenAra[3]);	useShdr = tokenAra[4];}catch (Exception e) {}	 
			    	scene.buildSierpinski(objName, scale, depth, useShdr != "No" );
			    	break;}
			    
			    //instance creation
			    case "named_object" : {
			    	//named_object <name>		    	
			    	String objName = tokenAra[1];
			    	scene.setObjectAsNamedObject(objName);			    	
			    	break;}
			    case "instance" : {
			    	//instance <name> <use current shader>
			    	String objName = tokenAra[1];
			    	String useCurShdr = "";				//whether to use the currently set shader for this instance, instead of the default shader for the named object
			    	if(tokenAra.length>2) {					useCurShdr = tokenAra[2];		    	}
			    	scene.addInstance(objName, useCurShdr != "");		    	
			    	break;}	
			    
			    /////////
			    // Textures
			    case "image_texture" :
			    case "texture" : {//load texture to be used for subsequent objects.  will be overridden by a surface command
			    	//First token could either be the side (top or bottom) or the textureName (in which case it is assumed to be the top texture)
			    	String side = tokenAra[1];
			    	String textureName = tokenAra[1];			    	
			    	if (side.toLowerCase().equals("bottom")){
			    		textureName = tokenAra[2]; 
			    		//if specified as bottom, assume bottom texture
			    		scene.currTextureBottom = ((my_procApplet)pa).loadImage(textureDir+textureName);
			    		scene.setHasGlblTxtrdBtm(true);
			    		win.getMsgObj().dispInfoMessage("myRTFileReader", "parseStringArray", "Bottom surface texture loaded");      }
			    	else {
			    		//if not specified then assume texture goes on top and texture name is specified in first token
			    		if (side.toLowerCase().equals("top")){  		  	textureName = tokenAra[2];   }
			    		scene.currTextureTop = ((my_procApplet)pa).loadImage(textureDir+textureName);
			    		scene.setHasGlblTxtrdTop(true);
			    		win.getMsgObj().dispInfoMessage("myRTFileReader", "parseStringArray", "Top surface texture loaded");
			    	} 
			    	scene.txtrType = 1;		//texture type is set to image/none
			    	break;}			    
			    
			    //procedural textures
			    //noise as txtr -> token[1] is scale value	
			    case "noise" : {    	scene.setNoise(Double.parseDouble(tokenAra[1]), tokenAra);    	break; }			    
			    //set colors used by procedural texture : <noise color tag> <color r g b> <-specify multiple times for each color
			    //need to put color commands/list after proc txtr type command in cli file, before object
			    case "noise_color" :{			    	scene.setTxtrColor(tokenAra);				    	break;   }			    
				//the following may just have <typ> or may have up to color scale, or may have all values - use defaults for all values not specified
				//<typ> <noise scale> <numOctaves> <turbMult> <pdMult x y z> <multByPI 1/0 1/0 1/0> <useFwdTransform 0/1> 
			    //		<rndomize colors colorScale - if present then true> <color mult> <num overlays - if present, otherwise 1>
			    case "marble" :
			    case "marble2" :
				case "stone" : 
			    case "wood" : 
			    case "wood2" : {   	scene.setTexture(tokenAra);   	break;    }			 
			    
			    //polygons		
			    //begin shape - defaults to triangle
			    case "begin" : {
			    	try {	vertType = tokenAra[1];   	} catch (Exception e) {	}      //triangle is default;  quad will be specified
			     	myVertCount = 0;
			      	if (vertType.equals("quad")){	myPoly = new myQuad(scene);} 
			      	else {				     		myPoly = new myTriangle(scene);}
			    	break;}
			    //texture_coord u v
			    case "texture_coord" : {
				//specifies the texture coordinate for the next vertex that is to be created for a triangle. - each "texture_coord" command will come before the corresponding "vertex" command in the .cli files
			    	((Base_PlanarObject) myPoly).setTxtrCoord(Double.parseDouble(tokenAra[1]),Double.parseDouble(tokenAra[2]), myVertCount); 
			    	break;}
			    case "vertex" : {
			    	((Base_PlanarObject) myPoly).setVert(Double.parseDouble(tokenAra[1]),Double.parseDouble(tokenAra[2]),Double.parseDouble(tokenAra[3]), myVertCount);
			    	myVertCount++;
			    	break;}
			    case "end" : {//end shape - reset vars to add new triangle, finalize triangle, add to arraylist of sceneobjects			    	
		    		((Base_PlanarObject) myPoly).finalizePoly();
		    		myPoly.shdr = scene.getCurShader();				//set shader for object <- important! this is where any proc textures used are specified
		    		scene.addObjectToScene(myPoly);
		    		vertType = "triangle";    //reset default object type as triangle
			    	myVertCount = 0;      
			    	break;}
			    			    
			    //prims - build in scene code
			    case "box" : 	    
			    case "plane" : 
			    case "cyl" :		   
			    case "cylinder" : 
			    case "hollow_cylinder" :
			    case "sphere" : 	    
			    case "moving_sphere":   
			    case "sphereIn" : 
			    case "ellipsoid" : {	readPrimData(scene, tokenAra);   	break;}
//				hollow_cylinder radius x z ymin ymax
//				Create a hollow cylinder that has its axis parallel to the y-axis. The hollow cylinder should not have end caps. 
			    
			    //matrix transformations of objects
			    case "push" : {  	scene.gtPushMatrix(); break;	    }//matrix push
			    case "pop" : {   	scene.gtPopMatrix();break;		    }//matrix pop
			    case "rotate" : {//needs to be in degrees as per assignment
			    	//builds a rotation matrix and adds to current transformation matrix on stack - angle, x,y,z
			    	scene.gtRotate(Double.parseDouble(tokenAra[1]), Double.parseDouble(tokenAra[2]), Double.parseDouble(tokenAra[3]), Double.parseDouble(tokenAra[4]));
			    	break;}
			    case "scale" : {
			    	//builds a scale matrix and adds to current transformation matrix on stack - sx,sy,sz
			    	scene.gtScale(Double.parseDouble(tokenAra[1]), Double.parseDouble(tokenAra[2]), Double.parseDouble(tokenAra[3]));
			    	break;}
			    case "translate" : {
			    	//builds a translate matrix and adds to current transformation matrix on stack - tx,ty,tz
			    	scene.gtTranslate(Double.parseDouble(tokenAra[1]), Double.parseDouble(tokenAra[2]), Double.parseDouble(tokenAra[3]));
			    	break;}	
			    			    
			    default : {
			    	win.getMsgObj().dispErrorMessage("myRTFileReader", "parseStringArray", "When reading "+fileName+" unknown command encountered : '"+tokenAra[0]+"' on line : ["+i+"] : " + fileStrings[i]);
			    	_printTokenAra(fileStrings[i], i, tokenAra);
			    	System.exit(0);
			    	break;}
			}//switch on object type from file	
		}//for each token string
		if((isMainFileScene) && !finalized){finalizeScene(loadedScenes,fileName, scene);	}	//puts finished scene -- only put in if _scene is null, meaning this is root file, not 2ndary read file
		return scene;
	}//parseStringArray
	
	/**
	 * Print contents of file and parsed results for each line
	 * @param fileStrings
	 */
	private void _debugFileStrings(String[] fileStrings) {
		for (int i=0; i<fileStrings.length; ++i) { 
			//Skip comments and empty lines
			if((fileStrings[i].startsWith("#")) || (fileStrings[i].strip().length() == 0)) {
				win.getMsgObj().dispInfoMessage("myRTFileReader", "_debugFileStrings", "\tLine " + i + ": fileString=`" + fileStrings[i] + "` : Comment or empty");
				continue;
			}
			String[] tokenAra = fileStrings[i].strip().split("\\s+"); // Get a line and parse tokens.
			//debug
			_printTokenAra(fileStrings[i], i, tokenAra);
		}
	}
	
	private void _printTokenAra(String fileString, int i, String[] token) {
		String tmpToken = "";
		for(String tok : token) {tmpToken += "|"+tok;}
		tmpToken += "|";
		win.getMsgObj().dispInfoMessage("myRTFileReader", "_printTokenAra", "\tLine " + i + ": fileString=`" + fileString + "` | tokens : "+ token.length + " : "+tmpToken);
	}//printTokenAra
	
	
	/**
	 * read in prim data from txt file and create object
	 * @param tokens
	 */
	private void readPrimData(Base_Scene scn, String[] tokens){
		Base_SceneObject tmp = null;
		switch(tokens[0]){
		    case "box" : {//box xmin ymin zmin xmax ymax zmax :
		    	double[] xAra = MyMathUtils.minAndMax(new double[]{Double.parseDouble(tokens[1]),Double.parseDouble(tokens[4])}),
		    			yAra = MyMathUtils.minAndMax(new double[]{Double.parseDouble(tokens[2]),Double.parseDouble(tokens[5])}),
		    			zAra = MyMathUtils.minAndMax(new double[]{Double.parseDouble(tokens[3]),Double.parseDouble(tokens[6])});
		    	double ctrX = (xAra[0] + xAra[1])*.5,
		    		ctrY = (yAra[0] + yAra[1])*.5,
		    		ctrZ = (zAra[0] + zAra[1])*.5;
		    	//putting box as a rendered bbox to minimize size of pure bboxes - rendered bbox is a bbox + shdr ref + some shdr-related functions and vars.
		    	tmp = new myRndrdBox(scn,ctrX, ctrY, ctrZ, new myVector(xAra[0], yAra[0], zAra[0]),	new myVector(xAra[1], yAra[1], zAra[1]));
		    	break;}			    
		    case "plane" : {			//infinite plane shape tkns 1-4 are A-D for plane equation (Ax + By + Cz + D = 0), token[5] is scale value to use for points (if not present default to 1)
		    	double scaleVal =  (tokens.length > 5) ? Double.parseDouble(tokens[5]) : 1.0;
		    	tmp = new myPlane(scn, Double.parseDouble(tokens[1]), Double.parseDouble(tokens[2]),Double.parseDouble(tokens[3]),Double.parseDouble(tokens[4]), scaleVal);
		    	break;}
		    case "cyl" : {//old cylinder code : cyl radius height center_x center_y center_z (orient_x orient_y orient_z : orientation is optional)
		    	double rad = Double.parseDouble(tokens[1]), hght = Double.parseDouble(tokens[2]);
		    	double xC = Double.parseDouble(tokens[3]), yC = Double.parseDouble(tokens[4]),zC = Double.parseDouble(tokens[5]);
		    	//check if orientation vector is given. Assume up (y==1)
		    	double xO = 0, yO = 1, zO = 0;
		    	if(tokens.length > 6) {
		    		xO = Double.parseDouble(tokens[6]);yO = Double.parseDouble(tokens[7]);zO = Double.parseDouble(tokens[8]);
		    	}
		    	tmp = new myCylinder(scn, rad,hght,xC,yC,zC,xO,yO,zO);
		    	break;}			   
		    case "cylinder" : { //alternative cylinder : cylinder radius x z ymin ymax  aligned on y axis
		    	double rad = Double.parseDouble(tokens[1]), xC = Double.parseDouble(tokens[2]), zC = Double.parseDouble(tokens[3]);
		    	double yMin  = Double.parseDouble(tokens[4]), yMax = Double.parseDouble(tokens[5]);
		    	double hght = yMax - yMin;
		    	//no given orientation vector. Assume up (y==1)
		    	double xO = 0,yO = 1,zO = 0;
		    	tmp = new myCylinder(scn, rad,hght,xC,yMin,zC,xO,yO,zO);
		    	break;}			    

		    case "hollow_cylinder" : {//hollow_cylinder radius x z ymin ymax aligned on y axis
		    	double rad = Double.parseDouble(tokens[1]), xC = Double.parseDouble(tokens[2]), zC = Double.parseDouble(tokens[3]);
		    	double yMin  = Double.parseDouble(tokens[4]), yMax = Double.parseDouble(tokens[5]);
		    	double hght = yMax - yMin;
		    	//no given orientation vector. Assume up (y==1)
		    	double xO = 0,yO = 1,zO = 0;
		    	tmp = new myHollow_Cylinder(scn, rad,hght,xC,yMin,zC,xO,yO,zO);		    	
		    	break;}
		    case "sphere" : {
		    	//create sphere
		    	tmp = new mySphere(scn, Double.parseDouble(tokens[1]),Double.parseDouble(tokens[2]),Double.parseDouble(tokens[3]),Double.parseDouble(tokens[4]),true);
		    	break;}			    
		    case "moving_sphere": {//moving_sphere radius x1 y1 z1 x2 y2 z2
		    	tmp = new myMovingSphere(scn, 
					Double.parseDouble(tokens[1]),Double.parseDouble(tokens[2]),Double.parseDouble(tokens[3]),Double.parseDouble(tokens[4]), 
					Double.parseDouble(tokens[5]),Double.parseDouble(tokens[6]),Double.parseDouble(tokens[7]), true);
		    	break;}			    
		    case "sphereIn" : {
		    	//create sphere with internal reflections - normals point in
		    	tmp = new mySphere(scn, Double.parseDouble(tokens[1]),Double.parseDouble(tokens[2]),Double.parseDouble(tokens[3]),Double.parseDouble(tokens[4]),true);
		    	tmp.setIsInverted(true);
		    	break;}	
		    case "ellipsoid" : {//create elliptical sphere with 3 radii elements in each of 3 card directions			    	
		    	tmp = new mySphere(scn, 
		    			Double.parseDouble(tokens[1]),Double.parseDouble(tokens[2]),Double.parseDouble(tokens[3]),
		    			Double.parseDouble(tokens[4]),Double.parseDouble(tokens[5]),Double.parseDouble(tokens[6]));
		    	break;}
		    default :{
		    	win.getMsgObj().dispErrorMessage("myRTFileReader", "readPrimData", "Object type not handled : "+ tokens[0]);
		    	return;
		    }
		}//switch
		//set shader and texture
		tmp.shdr = scn.getCurShader();
		//add object to scene
		scn.addObjectToScene(tmp);		
	}//readPrimData
	
	
	/**
	 * build a color value from a string array read in from a cli file.  stIdx is position in array where first color resides
	 * @param token
	 * @param stIdx
	 * @return
	 */
	private myRTColor readColor(String[] token, int stIdx){return new myRTColor(Double.parseDouble(token[stIdx]),Double.parseDouble(token[stIdx+1]),Double.parseDouble(token[stIdx+2]));}
	
	public void finalizeScene(TreeMap<String, Base_Scene> loadedScenes, String fileName,Base_Scene scene){
		//finalize scene
		scene.srcFileNames.push(fileName);
		loadedScenes.put(fileName, scene);
	}
	
	/**
	 * handle surface and shiny commands (very similar in layout but slight differences - shiny will use "simple" version of transmittance, not TIR
	 * @param scene
	 * @param token
	 * @param useSimple
	 */
	public void setSurfaceShiny(Base_Scene scene, String[] token, boolean useSimple){
		myRTColor cDiff = readColor(token,1);
		myRTColor cAmb = readColor(token,4);
		myRTColor cSpec = readColor(token,7);
		double phongExp = Double.parseDouble(token[10]);
		double kRefl = Double.parseDouble(token[11]);
		double kTrans =0, rfrIdx=0;
		myRTColor permClr = new myRTColor();
		scene.setHasGlblTxtrdTop(false);
		scene.setHasGlblTxtrdBtm(false);
		//See if surface description extends to transmission/refraction
		if (token.length > 12) {
			kTrans = Double.parseDouble(token[12]);  
			if (token.length > 13) {
				//perm
				rfrIdx = Double.parseDouble(token[13]);	
				//color of permiability/translucency - if has any will have all 3 values
				if (token.length > 14) {
					permClr.set(Double.parseDouble(token[14]),Double.parseDouble(token[15]),Double.parseDouble(token[16]));
				} else {
					permClr.set(rfrIdx,rfrIdx,rfrIdx);					
				}
			}
		}
		scene.setSurface(cDiff,cAmb,cSpec,phongExp,kRefl,kTrans,rfrIdx,permClr);

		if((useSimple) && ((kTrans > 0) || (rfrIdx > 0))){scene.setHasSimpleRefr(true);}			
	}//setSurfaceShiny

}//class myRTFileReader 


