//! Defines a structure for fast triangle mesh generation

use std::convert::Into;
use std::iter::FromIterator;
use std::collections::HashMap;
use crate::mesh::Mesh;
use crate::mesh::math::*;


type Vec2 = Vector2<f64>;

/// Structure intended to hold the points and faces during the generation
///
/// The generation functions are all methods of this class, and can be called independently of what is already in the shape, the common structure is only there to economize the memory reallocations.
/// The methods that are generating points and faces only from arguments (and so independently of what is already in the Shape) are using the name of the operation they represents (eg: `triangulation`), whereas the methods that use existing points or faces to operate on are called using verbs (eg: `triangulate`)
#[derive(Debug, Default)]
pub struct Shape {
	/// buffer of points
	pub points: Vec<Vec3>,		
	/// triangle indices in the buffer of points
	pub faces: Vec<[u32; 3]>,	
}


impl Into<Mesh> for Shape {
	fn into(self) -> Mesh {
		// TODO ideally we would do this:
		// Mesh::new(&self.faces, &self.points)
		// but the current Mesh new function is asking an owned storage, with splited vector and face components
		let mut faces = Vec::with_capacity(self.faces.len()*3);
		for face in self.faces {
			faces.extend_from_slice(&face);
		}
		let mut points = Vec::with_capacity(self.points.len()*3);
		for point in self.points {
			points.extend_from_slice((point.as_ref() as &[f64; 3]).as_ref());
		}
		Mesh::new(&faces, &points)
	}
}

impl Shape {
	/// create with empty buffers.
	pub fn new() -> Self {
		Shape::default()
	}
	
	/// return the extreme coordinates of the shape content.
	///
	pub fn boundingbox(&self) -> Option<[Vec3; 2]> {
		if self.points.is_empty()
			{ None }
		else {
			let mut extremes = [self.points[0], self.points[0]];
			for p in self.points.iter() {
				if      p.x < extremes[0].x		{ extremes[0].x = p.x; }
				else if p.y > extremes[1].x		{ extremes[1].x = p.x; }
				if      p.y < extremes[0].y		{ extremes[0].y = p.y; }
				else if p.y > extremes[1].y		{ extremes[1].y = p.y; }
				if      p.z < extremes[0].z		{ extremes[0].z = p.z; }
				else if p.z > extremes[1].z		{ extremes[1].z = p.z; }
			}
			Some(extremes)
		}
	}
	
	/// insert an other shape, simply including its points and faces.
	///
	pub fn merge(&mut self, other: &Self) -> &mut Self {
		let offset = self.points.len() as u32;
		self.points.extend_from_slice(&other.points);
		self.faces.extend(other.faces.iter().map(|f| [f[0]+offset, f[1]+offset, f[2]+offset]));
		self
	}
	
	/// Merge points that are too close.
	///
	/// if specified, `distance` is the approximate radius in which the points will me considered to be merged.
	/// # WANING: `distance` is approximative, points can be merged if they are up to `2*distance`
	///
	pub fn merge_doubles(&mut self, distance: Option<f64>) -> &mut Self {
		if self.points.is_empty() { return self; }
		
		// get the distance step
		let step = match distance {
			Some(distance)	=> distance,
			None => {
				let b = self.boundingbox().unwrap();
				(b[1]-b[0]).magnitude() * 1e-6
			}
		};
		
		// hashtable of redirections for merge
		let mut merges = HashMap::<u32, u32>::new();
		
		// generate a hash table to find the points quickly (this is truncated coordinates)
		let mut placements = HashMap::with_capacity(self.points.len());
		// handle truncation errors on the result of place
		let derives = [[0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,1,0], [0,1,1], [1,0,1], [1,1,1]];
		
		let place = |pt: &Vec3|	{ [(pt.x/step) as i64, (pt.y/step) as i64, (pt.z/step) as i64] };
		for (i,point) in self.points.iter().enumerate() {
			let placement = place(point);
			let mut pair_found = false;
			
			for derive in derives.iter() {
				let key = [placement[0]+derive[0], placement[1]+derive[1], placement[2]+derive[2]];
				if placements.contains_key(&key) {
					merges.insert(i as u32, placements[&key]);
					pair_found = true;
					break;
				}
			}
			if ! pair_found		{ placements.insert(placement, i as u32); }
		}
		
		// simplify multiple redirections
		for (start,mut target) in merges.iter() {
			while merges.contains_key(target) {
				target = &merges[target];
				assert_ne!(target, start);
			}
		}
		
		self.merge_points(&merges);
		self
	}
	
	/// merge some points.
	///
	/// the argument must be like `{point_to_merge:  point_to_merge_to}`.
	///
	pub fn merge_points(&mut self, merges: &HashMap<u32, u32>) -> &mut Self {
		let mut reindex = Vec::with_capacity(self.points.len());
		{
			let mut j = 0;
			for i in 0 .. self.points.len() {
				if ! merges.contains_key(&(i as u32)) 	{
					self.points[j] = self.points[i];
					j += 1;
				}
				reindex.push(j as u32);
			}
			self.points.truncate(j);
		}
		{
			let mut j = 0;
			for i in 0 .. self.faces.len() {
				let old = self.faces[i];
				self.faces[j] = [reindex[old[0] as usize], reindex[old[1] as usize], reindex[old[2] as usize]];
				let face = self.faces[j];
				if face[0] != face[1] || face[1] != face[2] || face[0] != face[2] {
					j += 1;
				}	
			}
			self.faces.truncate(j);
		}
		self
	}
	
	pub fn remove_unused(&mut self) -> &mut Self {
		let mut usage = vec![false; self.points.len()];
		for face in self.faces.iter() {
			for pt in face.iter() 		{ usage[*pt as usize] = true; }
		}
		
		let mut j = 0;
		for i in 0 .. self.points.len() {
			if usage[i]	{
				self.points[j] = self.points[i];
				j += 1;
			}
		}
		self.points.truncate(j);
		
		self
	}
	
	
	/// insert a surface defined by its outline.
	///
	/// The ouline is considered closed (ie. the last point is connected to the first).
	pub fn triangulation(&mut self, outline: &[Vec3]) -> &mut Self {
		let startindex = self.points.len();
		self.points.extend_from_slice(outline);
		self.triangulate((startindex as u32 .. (startindex+outline.len()) as u32).collect());
		self
	}
	
	/// create the triangles for a surface defined by its `outline`. 
	///
	/// The ouline is considered closed (ie. the last point is connected to the first)
	pub fn triangulate(&mut self, closed_outline: Vec<u32>) -> &mut Self {
		let mut outline = closed_outline;	// move it to mutable
		
		// fast-way for the trivial case
		if outline.len() == 3 {
			self.faces.push([outline[0], outline[1], outline[2]]);
		}
		
		// project all the points on the face approximate plane
		let (x,y,_) = base_from_points(&mut outline.iter().map(|i| self.points[*i as usize]));
		let mut pts = Vec::from_iter(	outline.iter()
											.map(|i| self.points[*i as usize])
											.map(|p| Vec2::new(p.dot(x), p.dot(y)))
											);
		
		'triangulation: while outline.len() > 2 {
		
			// find the thinest triangle along an edge from outline that doesn't contains other points
			let mut tinyest_index = 0;
			let mut tinyest_surf = 0.;
			'search: for i in 0 .. outline.len() {
				// adjacent points to i
				let i2 = if i==outline.len()-1 {0} else {i+1};
				let i1 = if i==0 {outline.len()-1} else {i-1};
				// adjacents segments to point i
				let c = pts[i];
				let a = pts[i2] - c;
				let b = pts[i1] - c;
				
				let surf = a.perp_dot(b);
				
				// if surface is negative, then pass (the surface is then outside the outline)
				if surf <= tinyest_surf { continue 'search; }
				
				let mat = Matrix2::from_cols(a, b).invert().unwrap();
				
				// check that there is not point of the ouline inside the triangle
				'interiors: for (j,&p) in pts.iter().enumerate() {
					if j==i1 || j==i || j==i2	{ continue 'interiors; }
					let params = mat * (p - c);
					let (u,v) = (params[0], params[1]);
					if	0. < u && u < 1. && 0. < v && v <= 1.-u && v != 1.
						{ continue 'search; }
				}
				tinyest_index = i;
				tinyest_surf = surf;
			}
			assert!(tinyest_surf > 0., format!("no more feasible triangles in {:?}", pts));
			
			// create the triangle and update the outline
			let i = tinyest_index;
			// adjacent points to i
			let i2 = if i==outline.len()-1 {0} else {i+1};
			let i1 = if i==0 {outline.len()-1} else {i-1};
			self.faces.push([outline[i2] as u32, outline[i] as u32, outline[i1] as u32]);
			outline.remove(i);
			pts.remove(i);
		}
		self
	}
	
	/// extrude a line into a lateral surface
	pub fn extrusion(&mut self, line: &[Vec3], displt: Vec3) -> &mut Self {
		self.points.reserve(line.len()*2);
		self.faces.reserve((line.len()-1)*2);
		
		let num = line.len() as u32;
		let istart = self.points.len() as u32;
		// get the points in
		for point in line {
			self.points.push(*point);
			self.points.push(point+displt);
		}
		// add the triangles
		for i in istart+1 .. istart+num {
			self.faces.push([i-1, i,   i+1]);
			self.faces.push([i-1, i+1, i+2]);
		}
		self
	}
	
	/// extrude a line into a lateral surface, proceeding by steps.
	///
	/// each point the result of `transform(step/segments, line_point)`
	///
 	pub fn extrans(&mut self, line: &[Vec3], segments: usize, transform: &dyn Fn(f64,Vec3) -> Vec3) -> &mut Self {
		self.points.reserve(line.len() * (segments+1));
		self.faces.reserve((line.len()-1) * segments * 2);
		
		let num = line.len() as u32;
		let mut istart = self.points.len() as u32;
		self.points.extend_from_slice(line);
		
		for segt in 1 ..= segments {
			let amount = segt as f64 / segments as f64;
			for pt in line.iter() {
				self.points.push(transform(amount, *pt));
			}
			for i in istart+1 .. istart+num {
				self.faces.push([i-1, i, i+num]);
				self.faces.push([i-1, i+num, i+num-1]);
			}
			istart += line.len() as u32;
		}
		self
 	}
 	
 	/// create a revolution surface from the line, around the axis.
 	///
 	pub fn revolution(&mut self, line: &[Vec3], segments: usize, origin: Vec3, axis: Vec3, angle: f64) -> &mut Self {
		self.extrans(line, segments, &|amount, pt| Quaternion::from_axis_angle(axis, Rad(angle*amount)).rotate_vector(pt-origin) + origin);
		self
 	}
 	
 	pub fn junction_shortest<J: JunctionMethod>(&mut self, line1: &[Vec3], line2: &[Vec3], method: J) -> &mut Self {
		let i1 = self.points.len() as u32;
		self.points.extend_from_slice(line1);
		let i2 = self.points.len() as u32;
		self.points.extend_from_slice(line2);
		let i3 = self.points.len() as u32;
		self.join_shortest(&Vec::from_iter(i1 .. i2), &Vec::from_iter(i2 .. i3), method);
		self
 	}
 	pub fn join_shortest<J: JunctionMethod>(&mut self, line1: &[u32], line2: &[u32], mut method: J) -> &mut Self {
		method.start(self, [0, 0]);
		let mut last1 = 0;
		let mut last2 = 0;
		while last1 < line1.len()-1 || last1 < line2.len()-1 {
			let next1 = if last1 < line1.len()-1  {last1+1} else {last1};
			let next2 = if last2 < line2.len()-1  {last2+1} else {last2};
			
			let dnext = (self.points[line1[next1] as usize] - self.points[line2[next2] as usize]).magnitude();
			let dlast1 = (self.points[line1[last1] as usize] - self.points[line2[next2] as usize]).magnitude();
			let dlast2 = (self.points[line1[next1] as usize] - self.points[line2[last2] as usize]).magnitude();
			if dnext < dlast1 && dnext < dlast2 {
				method.avance(self, [line1[next1], line2[next2]]);
				last1 = next1;
				last2 = next2;
			}
			else if dlast1 < dlast2 {
				method.avance(self, [line1[last1], line2[next2]]);
				last2 = next2;
			}
			else if dlast2 < dlast1 {
				method.avance(self, [line1[next1], line2[last2]]);
				last1 = next1;
			}
			else {
				panic!(format!("unable to determine the shortest between distances {} {} {}", dnext, dlast1, dlast2));
			}
		}
		self
 	}
 	 	
 	 	
 	pub fn junction_ratio<J: JunctionMethod>(&mut self, line1: &[Vec3], line2: &[Vec3], method: J) -> &mut Self {
		let i1 = self.points.len() as u32;
		self.points.extend_from_slice(line1);
		let i2 = self.points.len() as u32;
		self.points.extend_from_slice(line2);
		let i3 = self.points.len() as u32;
		self.join_ratio(&Vec::from_iter(i1 .. i2), &Vec::from_iter(i2 .. i3), method);
		self
 	}
 	pub fn join_ratio<J: JunctionMethod>(&mut self, line1: &[u32], line2: &[u32], mut method: J) -> &mut Self {
		method.start(self, [0, 0]);
		let length1: f64 = line1.windows(2).map(|w| (self.points[w[0] as usize] - self.points[w[1] as usize]).magnitude()).sum();
		let length2: f64 = line2.windows(2).map(|w| (self.points[w[0] as usize] - self.points[w[1] as usize]).magnitude()).sum();
		let mut last1 = 0;
		let mut last2 = 0;
		let mut absc1 = 0.;
		let mut absc2 = 0.;
		while last1 < line1.len()-1 || last1 < line2.len()-1 {
			let next1 = if last1 < line1.len()-1  {last1+1} else {last1};
			let next2 = if last2 < line2.len()-1  {last2+1} else {last2};
			
			let d1 = absc1 + (self.points[line1[last1] as usize] - self.points[line1[next1] as usize]).magnitude() / length1;
			let d2 = absc2 + (self.points[line2[last2] as usize] - self.points[line2[next2] as usize]).magnitude() / length2;
			if d1 < absc2 {
				method.avance(self, [line1[next1], line2[last2]]);
				last1 = next1;
				absc1 = d1;
			}
			else if d2 < absc1 {
				method.avance(self, [line1[last1], line2[next2]]);
				last2 = next2;
				absc2 = d2;
			}
			else {
				method.avance(self, [line1[next1], line2[next2]]);
				last1 = next1;
				last2 = next2;
				absc1 = d1;
				absc2 = d2;
			}
		}
		self
 	}
 	
 	fn create_tri(&mut self, pts: [u32; 3])	{ 
		self.faces.push(pts); 
	}
 	fn create_quad(&mut self, pts: [u32; 4]) { 
		self.faces.push([pts[0], pts[1], pts[2]]); 
		self.faces.push([pts[2], pts[3], pts[0]]); 
		// TODO: add quality checks (choose where to place the diagonal)
	}
 	
 	pub fn is_valid(&self) -> bool {
		// TODO use an error as return value
		let maxindex = self.points.len() as u32;
		for face in self.faces.iter() {
			for &i in face.iter() {
				if i > maxindex { return false; }
			}
		}
		true
 	}
}



/// method to create faces between 2 lines
pub trait JunctionMethod {
	/// called to create the first linking edge (pair bt two lines' points)
	/// the shape passed will then be always the same
	fn start(&mut self, shape: &mut Shape, pts: [u32; 2]);
	/// called to create faces between last linking edge and the given one
	fn avance(&mut self, shape: &mut Shape, pts: [u32; 2]);
}

/// very simple junction: create triangles between the two lines
pub struct JunctionSimple {
	lastedge: Option<[u32; 2]>,
}
impl JunctionSimple {
	fn new() -> Self {
		Self {lastedge: None}
	}
}
impl JunctionMethod for JunctionSimple {
	fn start(&mut self, _: &mut Shape, pts: [u32; 2])						{ self.lastedge = Some(pts); }
	fn avance(&mut self, shape: &mut Shape, pts: [u32; 2]) 	{
		let ptl = self.lastedge.unwrap();
		if      pts[0] == ptl[0]	{ shape.create_tri([pts[0], ptl[1], pts[1]]); }
		else if pts[1] == ptl[1]	{ shape.create_tri([pts[0], ptl[0], pts[1]]); }
		else						{ shape.create_quad([ptl[0], ptl[1], pts[1], pts[0]]); }
	}
}


/// Smooth junction: create a custom number of layers between the two lines.
/// The intermediate lines will be linear combinations of enclosing lines.
pub struct JunctionSmooth {
	segments: u32,
	lastedge: Option<([u32;2], u32)>,
}
impl JunctionSmooth {
	/// Create with a number of intermediate lines
	pub fn new(segments: u32) -> Self {
		Self {segments: segments, lastedge: None}
	}
	
	/// create the points for a linking edge, returing the start index of the new points
	fn create_edge(&self, shape: &mut Shape, pts: [u32;2]) -> u32 {
		let startindex = shape.points.len() as u32;
		let start = shape.points[pts[0] as usize];
		let end = shape.points[pts[1] as usize];
		for i in 0 .. self.segments {
			let f = ((i+1) as f64)/((self.segments+1) as f64);
			shape.points.push((f*start + (1.-f) * end));
		}
		startindex
	}
}
impl JunctionMethod for JunctionSmooth {
	fn start(&mut self, shape: &mut Shape, pts: [u32;2]) {
		self.lastedge = Some((pts, self.create_edge(shape, pts)));
	}

	fn avance(&mut self, shape: &mut Shape, pts: [u32; 2]) {
		let (ptl, edge1) = self.lastedge.unwrap();
		let edge2 = self.create_edge(shape, pts);
		self.lastedge = Some((pts, edge2));
		
		if self.segments > 0 {
			if ptl[0] == pts[0]		{ shape.create_tri([pts[0], edge1, edge2]); }
			else					{ shape.create_quad([ptl[0], edge1, edge2, pts[0]]); }
			for i in 0 .. self.segments-1 {
				shape.create_quad([edge1+i, edge1+i+1, edge2+i+1, edge2+i]);
			}
			if ptl[1] == pts[1]		{ shape.create_tri([edge2+self.segments-1, edge1+self.segments-1, ptl[1]]); }
			else 					{ shape.create_quad([edge2+self.segments-1, edge1+self.segments-1, ptl[1], pts[1]]); }
		}
		else {
			if      pts[0] == ptl[0]	{ shape.create_tri([pts[0], ptl[1], pts[1]]); }
			else if pts[1] == ptl[1]	{ shape.create_tri([pts[0], ptl[0], pts[1]]); }
			else						{ shape.create_quad([ptl[0], ptl[1], pts[1], pts[0]]); }
		}
	}
}


/*
pub struct JunctionSmoothDirs {
	segments: usize,
	lastedge: Option<([u32;2], u32)>,
	directions: HashMap<u32, (Vec3, f64)>,
}
*/

fn base_from_points(points: &mut dyn Iterator<Item=Vec3>) -> (Vec3, Vec3, Vec3) {
	let err = "not enough vectors to cet a base, needs 3 differents";
	let o = points.next().expect(err);
	let x = (points.next().expect(err) - o).normalize();
	let mut z = Vec3::zero();
	while z.magnitude() < 1e-4 {
		z = x.cross(points.next().expect(err) - o).normalize();
	}
	let y = z.cross(x);
	(x,y,z)
}


#[cfg(test)]
mod tests {
	use super::*;
	use std::f64::consts::PI;
	
	#[test]
	fn test_triangulation() 
	{
		let mut shape = Shape::new();
		shape.triangulation(&vec![
			Vec3::new(0., 0., 0.), 
			Vec3::new(1., 0., 0.), 
			Vec3::new(1., 1., 0.), 
			Vec3::new(2., 1., 0.),
			Vec3::new(0., 2., 0.),
			]);
		assert_eq!(shape.faces.len(), 3);
		assert!(shape.is_valid());
	}
	
	#[test]
	fn test_extrusion() {
		let mut shape = Shape::new();
		shape.extrusion(&vec![
			Vec3::new(0., 0., 0.), 
			Vec3::new(1., 0., 0.), 
			Vec3::new(1., 1., 0.), 
			Vec3::new(0., 1., 0.),
			], Vec3::new(0., 0.5, 1.));
		
		assert_eq!(shape.points.len(), 8);
		assert_eq!(shape.faces.len(), 6);	// there is no face for the non-existing edge (0, 3)
		assert!(shape.is_valid());
	}
	
	#[test]
	fn test_revolution() {
		let div = 32;
		let mut shape = Shape::new();
		shape.revolution(&vec![
			Vec3::new(1., 0., 0.),
			Vec3::new(2., 0., 1.),
			Vec3::new(1., 0., 2.),
			], div, Vec3::new(0., 0., 0.), Vec3::new(0., 0., 1.), PI);
		
		assert_eq!(shape.points.len(), 3*(div+1));
		assert_eq!(shape.faces.len(), 4*div);
		assert!(shape.is_valid());
	}
	
	#[test]
	fn test_merge_doubles() {
		let div = 32;
		let mut shape = Shape::new();
		shape.revolution(&vec![
			Vec3::new(1., 0., 0.),
			Vec3::new(2., 0., 1.),
			Vec3::new(1., 0., 2.),
			], div, Vec3::new(0., 0., 0.), Vec3::new(0., 0., 1.), 2.*PI);
		shape.merge_doubles(None);
		
		assert_eq!(shape.points.len(), 3*div);
		assert_eq!(shape.faces.len(), 4*div);
		assert!(shape.is_valid());
	}
	
	#[test]
	fn test_junction_simple() {
		let line1 = vec![
				Vec3::new(1., 0., 0.),
				Vec3::new(1., 1., 0.),
				Vec3::new(0., 1., 0.),
				];
		let line2 = vec![
				Vec3::new(1., 0., 1.),
				Vec3::new(1., 1., 1.),
				Vec3::new(0., 1., 1.),
				];
		{	// test shortest pairing
			let mut shape = Shape::new();
			shape.junction_shortest(&line1, &line2, JunctionSimple::new());
			assert_eq!(6, shape.points.len());
			assert_eq!(4, shape.faces.len());
			assert!(shape.is_valid());
		}
		{	// test ratio pairing
			let mut shape = Shape::new();
			shape.junction_ratio(&line1, &line2, JunctionSimple::new());
			assert_eq!(6, shape.points.len());
			assert_eq!(4, shape.faces.len());
			assert!(shape.is_valid());
		}
	}
	
	#[test]
	fn test_junction_smooth() {
		let line1 = vec![
				Vec3::new(1., 0., 0.),
				Vec3::new(2., 1., 0.),
				Vec3::new(0., 1., 0.),
				];
		let line2 = vec![
				Vec3::new(1., 0., 1.),
				Vec3::new(3., 1., 1.),
				Vec3::new(0., 2., 1.),
				];
		{	// test simple faces (reference test)
			let mut shape = Shape::new();
			shape.junction_shortest(&line1, &line2, JunctionSimple::new());
			assert_eq!(3*2, shape.points.len());
			assert_eq!(4*1, shape.faces.len());
		}
		{	// test smooth faces
			let mut shape = Shape::new();
			shape.junction_shortest(&line1, &line2, JunctionSmooth::new(10));
			assert_eq!(3*12, shape.points.len());
			assert_eq!(4*11, shape.faces.len());
			assert!(shape.is_valid());
		}
	}
}
