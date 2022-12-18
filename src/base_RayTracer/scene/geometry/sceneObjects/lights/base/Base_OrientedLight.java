package base_RayTracer.scene.geometry.sceneObjects.lights.base;

import base_Math_Objects.vectorObjs.doubles.myVector;
import base_RayTracer.scene.base.Base_Scene;
import base_RayTracer.scene.geometry.base.GeomObjType;

public abstract class Base_OrientedLight extends Base_Light{
	/**
	 * orientation of this light
	 */
	protected myVector orientation;

	public Base_OrientedLight(Base_Scene _scn, int _lightID, 
			double _r, double _g, double _b, 
			double _x, double _y, double _z, 
			double _dx, double _dy, double _dz, GeomObjType _type) {
		super(_scn,_lightID, _r, _g, _b, _x,_y,_z, _type);
		
	    orientation = new myVector(_dx,_dy,_dz); 
	    orientation._normalize();
	}
	
	public myVector getOrientation(double _t){	return orientation;	}
	/**
	 * send direction vector, finds multiplier for penumbra effect
	 * @param dir
	 * @param time
	 * @param innerThetRad
	 * @param outerThetRad
	 * @param radDiff
	 * @return
	 */
	protected double calcT_Mult(myVector dir, double time, double innerThetRad, double outerThetRad, double radDiff){
		double angle = Math.acos(-1*dir._dot(getOrientation(time)));			//intersection pt to light dir is neg light to intersection pt dir - want acos of this to get angle
		return getAngleProb(angle,innerThetRad,outerThetRad, radDiff);// (angle < innerThetRad ) ? 1 : (angle > outerThetRad ) ? 0 : (outerThetRad - angle)/radDiff;		
	}
	
}//class Base_OrientedLight
