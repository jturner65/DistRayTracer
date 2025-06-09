package base_RayTracer.scene.geometry.accelStruct.base;

import java.util.HashMap;
import java.util.Map;

/**
 * Enum defining the types of Acceleration Structure objects. Primarily BVH and Flat list now.
 * @author John Turner
 *
 */
public enum AccelStructType {	
	Unknown,
	FlatList,
	BVHTree,
	BVHLeftChild,
	BVHRightChild,
	BVHLeafList;
	
	private final String[] _typeExplanation = new String[]{
		"Unknown Acceleration Structure",
		"Flat List of Objects",
		"Bounding Volume Hierarchy Tree",
		"Bounding Volume Hierarchy Left Child",
		"Bounding Volume Hierarchy Right Child",
		"Bounding Volume Hierarchy Leaf List"			
	};
	private static final String[] _typeName = new String[]{
		"Unknown","Flat Object List","BVH Tree","BVH Left Child","BVH Right Child","BVH Leaf List"		
	};
	
	public static String[] getListOfTypes() {return _typeName;}
	private static Map<Integer, AccelStructType> map = new HashMap<Integer, AccelStructType>(); 
	static { for (AccelStructType enumV : AccelStructType.values()) { map.put(enumV.ordinal(), enumV);}}
	public int getOrdinal() {return ordinal();}
	public static AccelStructType getEnumByIndex(int idx){return map.get(idx);}
	public static AccelStructType getEnumFromValue(int idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
	public String getName() {return _typeName[ordinal()];}
	@Override
    public String toString() { return ""+_typeExplanation[ordinal()] + "("+ordinal()+")"; }	
    public String toStrBrf() { return ""+_typeExplanation[ordinal()]; }	
	
}//enum AccelStructType
