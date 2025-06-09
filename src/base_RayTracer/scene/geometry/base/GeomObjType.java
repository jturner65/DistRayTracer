package base_RayTracer.scene.geometry.base;

import java.util.HashMap;
import java.util.Map;

/**
 * Enum used to specify type of myGeomBase object
 * @author John Turner
 *
 */
public enum GeomObjType {
	None,BBox,RenderedBBox,Instance,PointLight,SpotLight,DiskLight,
	Triangle,Quad,Plane,Sphere,Cylinder,Hollow_Cylinder,Torus,
	AccelFlatList,AccelBVH;	

	private final String[] _typeExplanation = new String[]{
			"Non-object or unknown",
			"Bounding Box",
			"Renderable Bounding Box",
			"Instance of a Particular Object",
			"Point Light",
			"Spot Light",
			"Circular Area Light",
			"Triangle Planar Object",
			"Quadrilateral Planar Object",
			"Infinte Implicit Plane",
			"Implicit Sphere",
			"Implicit Capped Cylinder",
			"Implicit Cylindrical Tube",
			"Implicit Torus",
			"Flat List of Objects",
			"Bounding Volume Hierarchy Structure Holding Objects"};
	
	private static final String[] _typeName = new String[]{
			"None","Bounding Box","Rendered BBox", "Object Instance",
			"Point Light","Spot Light","Disk Light",
			"Triangle","Quad","Plane","Sphere","Capped Cylinder","Cylindrical Tube","Torus",
			"Flat List","BVH"};
	public static String[] getListOfTypes() {return _typeName;}
	private static Map<Integer, GeomObjType> map = new HashMap<Integer, GeomObjType>(); 
	static { for (GeomObjType enumV : GeomObjType.values()) { map.put(enumV.ordinal(), enumV);}}
	public int getOrdinal() {return ordinal();}
	public static GeomObjType getEnumByIndex(int idx){return map.get(idx);}
	public static GeomObjType getEnumFromValue(int idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
	public String getName() {return _typeName[ordinal()];}
	@Override
    public String toString() { return ""+_typeExplanation[ordinal()] + "("+ordinal()+")"; }	
    public String toStrBrf() { return ""+_typeExplanation[ordinal()]; }	

}//enum GeomObjType