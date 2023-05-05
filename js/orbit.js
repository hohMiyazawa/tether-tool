function roughFloatsAbsolute(float1,float2){
	return Math.abs(float1 - float2) < 0.00001;
}

function roughFloatsRelative(float1,float2){
	return roughFloatsAbsolute(float1/float2,1)
}

function roughFloatsAngles(float1,float2){
	const reduced1 = float1 % (2*Math.PI);
	const reduced2 = float2 % (2*Math.PI);
	return roughFloatsAbsolute(
		reduced1,
		reduced2
	) || roughFloatsAbsolute(
		Math.min(reduced1,2*Math.PI - reduced1) + Math.min(reduced2,2*Math.PI - reduced2),
		0
	)
}

function normaliseAngle(angle){
	let tmp = angle % (2*Math.PI);
	if(tmp < 0){
		return tmp + 2*Math.PI
	}
	return tmp
}

function bowlMin(func,low,high,iterations){
	let pivot1 = low;
	let pivot2 = high;
	for(let i=0;i<iterations;i++){
		let inter1 = pivot1 + (pivot2 - pivot1)/3;
		let inter2 = pivot1 + 2*(pivot2 - pivot1)/3;
		if(func(inter1) < func(inter2)){
			pivot2 = inter2
		}
		else{
			pivot1 = inter1
		}
	};
	return (pivot2 + pivot1)/2
}

const constants = {
	gravity: 6.6743e-11,
	boltzmann: 1.380649e-23,
	stefanBoltzmann: 5.670374419e-8,
	AU: 149597870700,
	c: 299792458,
	parsec: 648000*149597870700/Math.PI,
	lightyear: 299792458 * 86400,
}

function tisserand(orb1,orb2){
	return orb2.a/orb1.a + 2 * Math.cos(orb1.relativeInclination(orb2)) * Math.sqrt(orb1.a/orb2.a * (1 - Math.pow(orb1.e,2)))
}

let brok = function(num){
	if(num <= 0){
		if(num === 0){
			return [[0, 1, 0]]
		}
		return brok(-num).map(a => [-a[0],a[1],a[2]])
	}
	if(num >= 1){
		if(num === 1){
			return [[1, 1, 0]]
		}
		return brok(num % 1).map(a => [a[0] + a[1] * Math.trunc(num),a[1],a[2]])
	}
	let fracs = [];
	let lim1_a = 0;
	let lim1_b = 1;
	let lim2_a = 1;
	let lim2_b = 1;
	for(let i=0;i<1000;i++){
		if(lim1_a/lim1_b === num){
			fracs.push([lim1_a, lim1_b, 0]);
			return fracs
		}
		if(lim2_a/lim2_b === num){
			fracs.push([lim2_a, lim2_b, 0]);
			return fracs
		}
		let lowerDiff = num - lim1_a/lim1_b;
		let upperDiff = lim2_a/lim2_b - num;
		if(lowerDiff < upperDiff){
			if(
				!fracs.length
				|| fracs[fracs.length-1][0] !== lim1_a
				|| fracs[fracs.length-1][1] !== lim1_b
			){
				fracs.push([lim1_a, lim1_b, -lowerDiff])
			}
		}
		else{
			if(
				!fracs.length
				|| fracs[fracs.length-1][0] !== lim2_a
				|| fracs[fracs.length-1][1] !== lim2_b
			){
				fracs.push([lim2_a, lim2_b, upperDiff])
			}
		}
		let lim3_a = lim1_a + lim2_a;
		let lim3_b = lim1_b + lim2_b;
		if(lim3_a/lim3_b > num){
			lim2_a = lim3_a;
			lim2_b = lim3_b;
		}
		else{
			lim1_a = lim3_a;
			lim1_b = lim3_b;
		}
	}
	return fracs
}

const defaults = {
	geology: {
		k2: 0.2
	}
}

class Time{
	constructor(seconds){
		this.time = seconds
	}
	get seconds(){
		return this.time
	}
	get minutes(){
		return this.time/60
	}
	get hours(){
		return this.time/3600
	}
	get days(){
		return this.time/86400
	}
	get years(){
		return this.time/(365.25*86400)
	}
	valueOf(){
		return this.time
	}
}

class Vector{
	constructor(x,y,z){
		this.x = x;
		this.y = y;
		this.z = z;
	}
	get length(){
		return Math.hypot(this.x,this.y,this.z)
	}
	static crossProduct(vec1,vec2){
		return new Vector(
			vec1.y * vec2.z - vec1.z * vec2.y,
			vec1.z * vec2.x - vec1.x * vec2.z,
			vec1.x * vec2.y - vec1.y * vec2.x
		)
	}
	static dotProduct(vec1,vec2){
		return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z
	}
	dotProduct(vec){
		return Vector.dotProduct(this,vec)
	}
	static cos_theta(vec1,vec2){
		return vec1.dotProduct(vec2) / (vec1.length * vec2.length)
	}
}

class Transfer{
	constructor(properties){
		this.orbits = properties.orbits;
		this.burns = properties.burns;
		this.times = properties.times;
		this.name = properties.name || ""
	}
	get cost(){
		return this.burns.reduce((acc,val) => acc + val,0)
	}
	valueOf(){
		return this.cost
	}
}

let systems = new Map();

class System{
	constructor(properties){
		this.properties = properties
	}
	get orbit(){
		return this.properties.orbit
	}
	get atmosphere(){
		return this.properties.atmosphere
	}
	get geology(){
		return this.properties.geology
	}
	get GM(){
		return this.properties.GM || (this.properties.mass || Math.pow(this.radius,3)*4*Math.PI/3 * 2000) * constants.gravity
	}
	get mass(){
		return this.properties.mass || this.properties.GM / constants.gravity || Math.pow(this.radius,3)*4*Math.PI/3 * 2000
	}
	get radius(){
		let radius = this.properties.radius
			|| this.properties.diameter/2
			|| this.properties.radiusEquator
			|| this.properties.radiusPolar
			|| this.properties.diameterEquator/2
			|| this.properties.diameterPolar/2
		let radiusObject = {
			"polar": this.properties.radiusPolar || this.properties.diameterPolar/2 || radius,
			"equator": this.properties.radiusEquator || this.properties.diameterEquator/2 || radius,
			"valueOf": function(){return radius}
		};
		return radiusObject;
	}
	get diameter(){
		let diameter = this.properties.diameter
			|| this.properties.radius * 2
			|| this.properties.diameterPolar
			|| this.properties.diameterEquator
			|| this.properties.radiusPolar * 2
			|| this.properties.radiusEquator * 2
		let diameterObject = {
			"polar": this.properties.diameterPolar || this.properties.radiusPolar * 2 || diameter,
			"equator": this.properties.diameterEquator || this.properties.radiusEquator * 2 || diameter,
			"valueOf": function(){return diameter}
		};
		return diameterObject;
	}
	get volume(){
		if(this.properties.volume){
			return this.properties.volume
		}
		return this.radius.equator*this.radius.equator*this.radius.polar*4*Math.PI/3
	}
	get area(){
		if(this.properties.area){
			return this.properties.area
		}
		if(this.radius.equator === this.radius.polar){//sphere
			return 4*Math.PI*this.radius*this.radius;
		}
		else if(this.radius.equator < this.radius.polar){//prolate
			let eccentricity = Math.sqrt(1 - Math.pow(this.radius.equator,2)/Math.pow(this.radius.equator,2));
			return 2*Math.PI*Math.pow(this.radius.equator,2)*(1 + this.radius.polar*Math.asin(eccentricity)/(this.radius.equator*this.radius.polar))
		}
		else if(this.radius.equator > this.radius.polar){//oblate
			let eccentricity = Math.sqrt(1 - Math.pow(this.radius.polar,2)/Math.pow(this.radius.equator,2));
			return 2*Math.PI*Math.pow(this.radius.equator,2) + Math.PI*Math.pow(this.radius.polar,2)*Math.log((1 + eccentricity)/(1 - eccentricity))/eccentricity;
		}
	}
	get density(){
		return this.properties.density
			|| this.mass / this.volume
	}
	get gravity(){
		return this.properties.gravity
			|| this.GM / Math.pow(this.radius,2)
	}
	get period(){
		return new Time(this.properties.period)
	}
	get surfaceEscape(){
		return Math.sqrt(2*this.GM/this.radius)
	}
	get surfaceVelocity(){
		return Math.PI*this.diameter/this.period
	}
	get tidalLock(){
		return this.properties.tidalLock === true
	}
	get satellites(){
		return (this.properties.satellites || []).map(e => systems.get(e) || e)
	}
	get stationaryOrbit(){
		return new Orbit({
			system: this,
			radius: Math.cbrt(this.GM/Math.pow(2*Math.PI/this.period,2))
		})
	}
	get color(){
		return this.properties.color || "#BEBEBE"
	}
	get L1(){
		if(!this.properties.L1){
			let pivot1 = this.radius.equator;
			let pivot2 = this.orbit.a/2;
			const span = pivot2 - pivot1;
			const w_square = Math.pow(2*Math.PI/this.orbit.period,2);
			const mu1 = this.orbit.system.GM;
			const mu2 = this.GM;
			const a = this.orbit.a;
			while(pivot2 - pivot1 > span/1e12){
				let centre = (pivot1 + pivot2)/2;
				let force = mu2/Math.pow(centre,2) + w_square*(mu1/(mu1 + mu2) * a - centre) - mu1/Math.pow(a - centre,2);
				if(force > 0){
					pivot1 = centre
				}
				else{
					pivot2 = centre
				}
			}
			this.properties.L1 = (pivot1 + pivot2)/2
		}
		return this.properties.L1
	}
	get L2(){
		if(!this.properties.L2){
			let pivot1 = this.radius.equator;
			let pivot2 = this.orbit.a;
			const span = pivot2 - pivot1;
			const w_square = Math.pow(2*Math.PI/this.orbit.period,2);
			const mu1 = this.orbit.system.GM;
			const mu2 = this.GM;
			const a = this.orbit.a;
			while(pivot2 - pivot1 > span/1e12){
				let centre = (pivot1 + pivot2)/2;
				let force = - mu2/Math.pow(centre,2) + w_square*(mu1/(mu1 + mu2) * a + centre) - mu1/Math.pow(a + centre,2);
				if(force < 0){
					pivot1 = centre
				}
				else{
					pivot2 = centre
				}
			}
			this.properties.L2 = (pivot1 + pivot2)/2
		}
		return this.properties.L2
	}
	get lowOrbit(){
		let radius = this.radius.equator + (this.properties.lowOrbitAltitude || this.radius.equator/30);
		if(this.atmosphere && this.atmosphere.height){
			radius = Math.max(radius,this.radius.equator + this.atmosphere.height)
		}
		return new Orbit({
			system: this,
			radius: radius
		})
	}
	turningAngle(periapsis,velocity){
		return -2*Math.asin(-1/(1 + periapsis*Math.pow(velocity,2)/this.GM))
	}
	gravity(r){
		return this.GM/(r*r)
	}
	toString(){
		return this.properties.name || "unnamed"
	}
}

class Atmosphere{
	constructor(properties){
		Object.keys(properties).forEach(key => {
			this[key] = properties[key]
		});
		this.layers.forEach((layer,index) => {
			if(!layer.hasOwnProperty("lapseRate")){
				if(this.layers.length === 1){
					layer.lapseRate = 0
				}
				else if((index + 1) === this.layers.length){
					layer.lapseRate = this.layers[index - 1].lapseRate
				}
				else{
					layer.lapseRate = (layer.temperature - this.layers[index + 1].temperature)/((this.layers[index + 1].base - layer.base)/1000)
				}
			}
		})
	}
	pressure(alt){
	}
	temperature(alt){
		let layer = this.layers.filter(e => e.base < alt).pop();
		return layer.temperature - layer.lapseRate * (alt - layer.base)/1000
	}
	pressure(alt){
		let layer = this.layers.filter(e => e.base < alt).pop();
		return layer.pressure
	}
}

class Geology{
	constructor(properties,parent){
		this.properties = properties;
		this.parent = parent
	}
	get k2(){
		return this.properties.k2 || defaults.geology.k2
	}
	get tidalHeating_loveless(){
		return 21/2
			* constants.gravity
			* Math.pow(this.parent.orbit.system.mass,2)
			* Math.pow(this.parent.radius,5)
			* this.parent.orbit.meanMotion
			* Math.pow(this.parent.orbit.eccentricity,2)
			/ Math.pow(this.parent.orbit.a,6)
	}
	get tidalHeating(){
		return 21/2
			* constants.gravity
			* this.k2
			* Math.pow(this.parent.orbit.system.mass,2)
			* Math.pow(this.parent.radius,5)
			* this.parent.orbit.meanMotion
			* Math.pow(this.parent.orbit.eccentricity,2)
			/ Math.pow(this.parent.orbit.a,6)
	}
}
//new Orbit() <-- unit system with radius 0
//new Orbit(properties)
//new Orbit(radius)
//new Orbit(periapsis[,apoapsis[,longitude[,inclination[,argument[,anomaly]]]]])
class Orbit{
	highPrecision = {}
	constructor(propertiesOrPeriapsis,apoapsis,longitude,inclination,argument,anomaly){
		if(arguments.length === 0){
			this.system = systems.get("UNIT");
			this.periapsis = 1;
			this.apoapsis = 1;
			this.longitude = 0;
			this.inclination = 0;
			this.argument = 0;
			this.anomaly = 0;
		}
		else if(typeof propertiesOrPeriapsis === "number"){
			this.system = systems.get("UNIT");
			this.periapsis = propertiesOrPeriapsis;
			this.apoapsis = apoapsis || this.periapsis;
			this.longitude = longitude || 0;
			this.inclination = inclination || 0;
			this.argument = argument || 0;
			this.anomaly = anomaly || 0;
		}
		else{
			if(arguments.length > 1){
				throw "too many arguments";
			}
			let properties = propertiesOrPeriapsis;
			if(properties.system){
				if(typeof properties.system === "string"){
					this.system = systems.get(properties.system);
				}
				else{
					this.system = properties.system;
				}
			}
			else{
				this.system = systems.get("UNIT");
			}
			this.periapsis = properties.periapsis || properties.radius || properties.apoapsis || properties.semiMajor*(1 - properties.eccentricity) || 1;
			this.apoapsis = properties.apoapsis || properties.radius || properties.periapsis || properties.semiMajor*(1 + properties.eccentricity) || 1;
			this.longitude = properties.longitude || 0;
			this.inclination = properties.inclination || 0;
			this.argument = properties.argument || 0;
			this.anomaly = properties.anomaly || 0;
			if(properties.hasOwnProperty("semiMajor")){
				this.highPrecision.semiMajor = properties.semiMajor
			}
			if(properties.hasOwnProperty("period")){
				this.highPrecision.period = properties.period
			}
			if(properties.hasOwnProperty("eccentricity")){
				this.highPrecision.eccentricity = properties.eccentricity
			}
		}
		if(arguments.length > 6){
			throw "too many arguments";
		}
		if(Math.abs(this.apoapsis) < this.periapsis){
			throw "invalid orbit";
		}
	}
	get A(){
		return this.apoapsis
	}
	get P(){
		return this.periapsis
	}
	get semiMajor(){
		if(this.highPrecision.hasOwnProperty("semiMajor")){
			return this.highPrecision.semiMajor
		}
		return (this.A + this.P)/2
	}
	get a(){
		return this.semiMajor
	}
	radius(anomaly){
		return -2*this.apoapsis*this.periapsis/(
			this.apoapsis * (- Math.cos(anomaly) - 1)
			+ this.periapsis * (Math.cos(anomaly) - 1)
		)
	}
	velocityFromAnomaly(anomaly){
		let distance = this.radius(anomaly);
		let velocityObject = this.velocity(distance)
		if(normaliseAngle(anomaly) < Math.PI){
			return velocityObject
		}
		return {
			"tangential": velocityObject.tangential,
			"radial": -velocityObject.radial,
			"valueOf": function(){return velocityObject.valueOf()}
		}
	}
	get altitude(){
		return this.semiMajor - this.system.radius
	}
	get eccentricity(){
		if(this.highPrecision.hasOwnProperty("eccentricity")){
			return this.highPrecision.eccentricity
		}
		return (this.A - this.P)/(this.A + this.P)
	}
	get e(){
		return this.eccentricity
	}
	get type(){
		if(this.e === 0){
			return "circular"
		}
		else if(this.e < 1){
			return "elliptic"
		}
		else if(this.e === 1){
			return "parabolic"
		}
		else{
			return "hyperbolic"
		}
	}
	squaredVelocity(r){
		return this.system.GM*(2/r - 1/this.a)
	}
	get angularMomentum(){
		return Math.sqrt(this.squaredVelocity(this.P)) * this.P
	}
	get meanMotion(){
		return 2*Math.PI/this.period
	}
	velocity(r){
		let velocity = Math.sqrt(this.squaredVelocity(r));
		let velocityObject = {
			"tangential": this.angularMomentum/r,
			"radial": Math.sqrt(this.squaredVelocity(r) - Math.pow(this.angularMomentum/r,2)),
			"valueOf": function(){return velocity}
		};
		if(r === this.P || r === this.A){
			velocityObject.radial = 0
		}
		return velocityObject;
	}
	v(r){
		return this.velocity(r)
	}
	get velP(){
		let info = this.v(this.P);
		info.radial = 0;
		return info;
	}
	get velA(){
		let info = this.v(this.A);
		info.radial = 0;
		return info;
	}
	get vela(){
		return this.v(this.a);
	}
	get E(){
		return (this.velP*this.velP - 2/this.P)/2
	}
	get period(){
		let period = 2*Math.PI*Math.sqrt(Math.pow(this.a,3)/this.system.GM);
		if(this.highPrecision.hasOwnProperty("period")){
			period = this.highPrecision.period
		}
		return new Time(period)
	}
	get T(){
		return this.period
	}
	circular(r){
		let velocity = Math.sqrt(this.system.GM/r);
		return {
			"tangential": Math.sqrt(this.system.GM/r),
			"radial": 0,
			"valueOf": function(){return velocity}
		};
	}
	escape(r){
		return Math.sqrt(2*this.system.GM/r)
	}
	get vinf(){
		return Math.sqrt(this.squaredVelocity(this.P) - this.escape(this.P)*this.escape(this.P))
	}
	get escapeCost(){
		const cost = Math.abs(this.escape(this.P) - this.velP);
		return {
			"orbit": new Orbit({
				system: this.system,
				periapsis: this.P,
				apoapsis: Infinity,
				eccentricity: 1,
				inclination: this.inclination,
				longitude: this.longitude,
				anomaly: 0
			}),
			"valueOf": function(){return cost}
		}
	}
	inner(r){
		let a = this.a;let P = this.P;let A = this.A;
		//let grunn = a - (A + P)*(r - P)/(A - P);
		let grunn = (A + P)*(r - P)/(A - P);
		//grunn / (A + P) = (r - P)(A - P)
		let area = Math.PI*a*a;
		let sector = area * Math.acos((a - grunn)/a)/(2*Math.PI);
		let height = Math.sqrt(a*a - Math.pow(grunn - a,2));
		let triangle = height*(a - P)/2;
		return {
			angle: Math.PI - Math.acos((grunn - P)/r),
			time: new Time(this.period*(sector - triangle)/area),
			velocity: this.velocity(r)
		};
	}
	static vectorsToOrbit(system,position,velocity){
		if(typeof system === "string"){
			system = systems.get(system)
		}
		if(Array.isArray(position)){
			position = {
				x: position[0],
				y: position[1],
				z: position[2]
			}
		}
		if(Array.isArray(velocity)){
			velocity = {
				x: velocity[0],
				y: velocity[1],
				z: velocity[2]
			}
		};
		let radius = Math.hypot(position.x,position.y,position.z);
		let velMagnitude = Math.hypot(velocity.x,velocity.y,velocity.z);
		let energy = Math.pow(velMagnitude,2)/2 - system.GM/radius;
		let a = - 1 /(velMagnitude * velMagnitude / system.GM - 2/radius);
		let velTan = velMagnitude * Math.sqrt(1 - Math.pow(Vector.cos_theta(new Vector(position.x,position.y,position.z), new Vector(velocity.x,velocity.y,velocity.z)),2));
/*
v^2 = mu * (2/r - 1/a)
v^2 / mu= 2/r - 1/a
v^2 / mu - 2/r = - 1/a
-(v^2 / mu - 2/r) = 1/a
a = -1 / (v^2 / mu - 2/r) 


vp^2 = mu * (2/P - 1/a)
(v*r/P)^2 = mu * (2/P - 1/a)

*/

		let P = (a*system.GM - Math.sqrt(a*a * system.GM*system.GM - a * system.GM * radius * radius * velTan * velTan))/system.GM;
		let orbit = new Orbit({system: system, periapsis: P, apoapsis: 2*a - P});
		return orbit;
	}
	get planeVector(){
		const mag = Math.sin(this.inclination);
		return new Vector(
			mag * Math.sin(this.longitude),
			- mag * Math.cos(this.longitude),
			Math.cos(this.inclination)
		)
	}
	static identical(orbit1,orbit2){//this can't be guaranteed by numerics, but it's useful in some cases
		return orbit1.P === orbit2.P
			&& orbit1.A === orbit2.A
			&& (orbit1.inclination % (2*Math.PI) === orbit2.inclination % (2*Math.PI))
			&& (
				(orbit1.longitude % (2*Math.PI) === orbit2.longitude % (2*Math.PI))
				|| orbit1.inclination === 0
			)
			&& (
				(orbit1.argument % (2*Math.PI) === orbit2.argument % (2*Math.PI))
				|| (orbit1.A === orbit1.P)
			)
	}
	identical(orbit){
		return Orbit.identical(this,orbit)
	}
	static roughlyEqual(orbit1,orbit2){
		return roughFloatsRelative(orbit1.P,orbit2.P)
			&& roughFloatsRelative(orbit1.A,orbit2.A)
			&& roughFloatsAngles(orbit1.inclination,orbit2.inclination)
			&& (
				roughFloatsAngles(orbit1.longitude,orbit2.longitude)
				|| roughFloatsAngles(orbit1.inclination,0)
			)
			&& (
				roughFloatsAngles(orbit1.argument,orbit2.argument)
				|| roughFloatsRelative(orbit1.A,orbit1.P)
			)
	}
	roughlyEqual(orbit){
		return Orbit.roughlyEqual(this,orbit)
	}
	static relativeInclination(orbit1,orbit2){
		return Math.acos(
			Vector.dotProduct(
				orbit1.planeVector,
				orbit2.planeVector
			)
		)
	}
	relativeInclination(orbit){
		return Orbit.relativeInclination(this,orbit)
	}
	relativeNodeLine(orbit){
		let cross = Vector.crossProduct(
			this.planeVector,
			orbit.planeVector
		);
		let nodeVector = new Vector(
			Math.cos(this.longitude),
			Math.sin(this.longitude),
			0
		);
		if(cross.length === 0){
			return 0
		}
		let angle = Math.acos(Vector.dotProduct(cross,nodeVector)/cross.length);
		return angle - this.argument
	}
	static singleBurn(orbit1,orbit2){
		if(orbit1.apoapsis < orbit2.periapsis || orbit2.apoapsis < orbit1.periapsis){
			return null
		}
		let relativeInclination = orbit1.relativeInclination(orbit2);
		if(relativeInclination === 0){
			//TODO
		}
		else{
			let anomaly1 = normaliseAngle(orbit1.relativeNodeLine(orbit2));
			let anomaly2 = normaliseAngle(orbit2.relativeNodeLine(orbit1));
			let distance1a = orbit1.radius(anomaly1);
			let distance1b = orbit2.radius(anomaly2 + Math.PI);
			let distance2a = orbit1.radius(anomaly1 + Math.PI);
			let distance2b = orbit2.radius(anomaly2);
			let best1 = null;
			let best2 = null;
			if(roughFloatsRelative(distance1a,distance1b)){
				let vel1 = orbit1.velocityFromAnomaly(anomaly1);
				let vel2 = orbit2.velocityFromAnomaly(anomaly2 + Math.PI);
				best1 = new Transfer({
					name: "single-burn",
					orbits: [],
					times: [0],
					burns: [Math.hypot(
						vel1.radial - vel2.radial,
						vel1.tangential - vel2.tangential * Math.cos(relativeInclination),
						Math.sin(relativeInclination)*vel2.tangential
					)]
				})
			}
			if(roughFloatsRelative(distance2a,distance2b)){
				let vel1 = orbit1.velocityFromAnomaly(anomaly1 + Math.PI);
				let vel2 = orbit2.velocityFromAnomaly(anomaly2);
				best2 = new Transfer({
					name: "single-burn",
					orbits: [],
					times: [0],
					burns: [Math.hypot(
						vel1.radial - vel2.radial,
						vel1.tangential - vel2.tangential * Math.cos(relativeInclination),
						Math.sin(relativeInclination)*vel2.tangential
					)]
				})
			}
			if(best2 !== null && best1 !== null && best2 < best1){
				return best2
			}
			return best1 || best2
		}
	}
	static biTouch(orbit1,orbit2){
		let relativeInclination = orbit1.relativeInclination(orbit2);
		if(relativeInclination !== 0){
			return null
		}
		let best;
		let bestIndex = 0;
		const iterations = 64;
		if(orbit1.periapsis === orbit1.apoapsis && orbit2.periapsis === orbit2.apoapsis){
			return Orbit.hohmann(orbit1,orbit2)
		}
		else if(orbit1.periapsis === orbit1.apoapsis){
			if(orbit1.periapsis === orbit2.apoapsis){
				return Orbit.singleBurn(orbit1,orbit2)
			}
		}
		else if(orbit2.periapsis === orbit2.apoapsis){
		}
		let getTransfer = function(anomaly){
			let distance1 = orbit1.radius(anomaly - orbit1.argument);
			let distance2 = orbit2.radius(anomaly + Math.PI - orbit2.argument);
			let vel1 = orbit1.velocityFromAnomaly(anomaly - orbit1.argument);
			let vel2 = orbit2.velocityFromAnomaly(anomaly + Math.PI - orbit2.argument);
			let transferOrbit = new Orbit({
				system: orbit1.system,
				inclination: orbit1.inclination,
				periapsis: Math.min(distance1,distance2),
				apoapsis: Math.max(distance1,distance2),
				longitude: orbit1.longitude,
				argument: (distance1 < distance2 ? anomaly : anomaly + Math.PI)
			})
			let trans = new Transfer({
				name: "bi-touch",
				orbits: [transferOrbit],
				times: [transferOrbit.period/2],
				burns: [
					Math.hypot(vel1.radial,transferOrbit.velocity(distance1) - vel1.tangential),
					Math.hypot(vel2.radial,transferOrbit.velocity(distance2) - vel2.tangential)
				]
			})
			return trans
		}
		for(let i=0;i<iterations;i++){
			let trans = getTransfer(2*Math.PI*i/iterations);
			if(!best || best > trans){
				best = trans;
				bestIndex = i;
			}
		}
		for(let i=0;i<iterations;i++){
			let trans = getTransfer(
				2*Math.PI*(
					(
						bestIndex-1 + 2*i/(iterations - 1)
					)/iterations
				)
			);
			if(!best || best > trans){
				best = trans;
			}
		}
		return best
	}
	biTouch(orbit){
		return Orbit.biTouch(this,orbit)
	}
	static planeSplit(orbit1,orbit2){
		let anomaly1 = orbit1.relativeNodeLine(orbit2);
		let anomaly2 = orbit2.relativeNodeLine(orbit1);
		let distance1 = orbit1.radius(anomaly1);
		let distance2 = orbit2.radius(anomaly2);
		let relativeInclination = orbit1.relativeInclination(orbit2);
		if(relativeInclination === 0){
			return Orbit.biTouch(orbit1,orbit2)
		}

		let intermediate1 = new Orbit({
			system: orbit1.system,
			periapsis: Math.min(distance1,distance2),
			apoapsis: Math.max(distance1,distance2)
		});

		let vel1 = orbit1.velocity(distance1);
		let vel2 = orbit2.velocity(distance2);

		let angle1Cost = angle => {
			let firstTransfer = Math.hypot(
				vel1.radial,
				Math.cos(angle) * intermediate1.velocity(distance1) - vel1.tangential,
				Math.sin(angle) * intermediate1.velocity(distance1)
			)
			let secondTransfer = Math.hypot(
				vel2.radial,
				Math.cos(relativeInclination - angle) * intermediate1.velocity(distance2) - vel2.tangential,
				Math.sin(relativeInclination - angle) * intermediate1.velocity(distance2)
			);
			return firstTransfer + secondTransfer
		}
		let gamma = bowlMin(
			angle1Cost,
			0,
			relativeInclination,
			20
		);
		let transferCost = angle1Cost(gamma);
		anomaly1 += Math.PI;
		anomaly2 += Math.PI;
		distance1 = orbit1.radius(anomaly1);
		distance2 = orbit2.radius(anomaly2);
		let intermediate2 = new Orbit({
			system: orbit1.system,
			periapsis: Math.min(distance1,distance2),
			apoapsis: Math.max(distance1,distance2)
		});
		vel1 = orbit1.velocity(distance1);
		vel2 = orbit2.velocity(distance2);
		let angle2Cost = angle => {
			let firstTransfer2 = Math.hypot(
				vel1.radial,
				Math.cos(angle) * intermediate2.velocity(distance1) - vel1.tangential,
				Math.sin(angle) * intermediate2.velocity(distance1)
			)
			let secondTransfer2 = Math.hypot(
				vel2.radial,
				Math.cos(relativeInclination - angle) * intermediate2.velocity(distance2) - vel2.tangential,
				Math.sin(relativeInclination - angle) * intermediate2.velocity(distance2)
			)
			return firstTransfer2 + secondTransfer2
		}
		let gamma2 = bowlMin(
			angle2Cost,
			0,
			relativeInclination,
			20
		);
		let transferCost2 = angle2Cost(gamma2);
		if(transferCost < transferCost2){
			return new Transfer({
				name: "plane-split",
				burns: [transferCost],
				orbits: [intermediate1]
			})
		}
		else{
			return new Transfer({
				name: "plane-split",
				burns: [transferCost2],
				orbits: [intermediate2]
			})
		}
	}
	planeSplit(orbit){
		return Orbit.planeSplit(this,orbit)
	}
	static hohmann(orbit1,orbit2){
		if((typeof orbit1 === "number") && (typeof orbit2 === "number")){
			if(orbit1 > orbit2){
				[orbit1,orbit2] = [orbit2,orbit1]
			};
			let transfer = new Orbit(orbit1,orbit2);
			return new Transfer({
				name: "hohmann",
				orbits: [transfer],
				times:  [transfer.period/2],
				burns:  [transfer.velP - transfer.circular(orbit1),transfer.circular(orbit2) - transfer.velA]
			})
		}
		else if((typeof orbit1 !== "number") && (typeof orbit2 !== "number")){
			if(orbit1.system.name !== orbit2.system.name){
				return null
			}
			if(!roughFloatsAbsolute(orbit1.e,0) || !roughFloatsAbsolute(orbit2.e,0)){
				return null
			};
			if(!roughFloatsAbsolute(orbit1.inclination,orbit2.inclination)){
				return null
			};
			if(!roughFloatsAbsolute(orbit1.longitude,orbit2.longitude)){
				return null
			};
			if(orbit1.a > orbit2.a){
				[orbit1,orbit2] = [orbit2,orbit1]
			};
			let transfer = new Orbit({
				system: orbit1.system,
				periapsis: orbit1.a,
				apoapsis: orbit2.a
			});
			return new Transfer({
				name: "hohmann",
				orbits: [transfer],
				times: [transfer.period/2],
				burns: [transfer.velP - transfer.circular(orbit1.a),transfer.circular(orbit2.a) - transfer.velA]
			})
		}
		else{
			if(typeof orbit1 === "number"){
				orbit1 = new Orbit({
					system: orbit2.system,
					radius: orbit1
				})
			}
			else{
				orbit2 = new Orbit({
					system: orbit1.system,
					radius: orbit2
				})
			}
			return Orbit.hohmann(orbit1,orbit2)
		}
	}
	hohmann(orbit){
		return Orbit.hohmann(this,orbit)
	}
	static biElliptic(orbit1,orbit2){
		return new Transfer({
			name: "bi-elliptic",
			orbits: [],
			burns: [orbit1.escapeCost,orbit2.escapeCost]
		})
	}
	biElliptic(orbit){
		return Orbit.biElliptic(this,orbit)
	}
	static nodeToNode(orbit1,orbit2){
		let relative = orbit1.relativeInclination(orbit2)
	}
	nodeToNode(orbit){
		return Orbit.nodeToNode(orbit)
	}
	static transfer(orbit1,orbit2){
		if(orbit1.system.properties.name !== orbit2.system.properties.name){
			if(orbit1.system.orbit.system.properties.name == orbit2.system.orbit.system.properties.name){
				let interplanetary = orbit1.system.orbit.planeSplit(orbit2.system.orbit);
				console.log(interplanetary);
			}
			else{
				console.warn("patched conics transfer not implemented");
				return null;
			}
		}
		else{
			let best = Orbit.biElliptic(orbit1,orbit2);
			let splitPlaneCost = orbit1.planeSplit(orbit2);
			if(splitPlaneCost < best){
				best = splitPlaneCost
			}
			return best
		}
	}
	transfer(orbit){
		return Orbit.transfer(this,orbit)
	}
	reachRadius(r){
		if(r >= this.P && (this.A < 0 || this.A >= r)){
			return {
				"orbit": new Orbit(this),
				"valueOf": function(){return 0}
			}
		}
		if(r > this.A){
			let raisedOrbit = new Orbit({
				system: this.system,
				periapsis: this.periapsis,
				apoapsis: r,
				inclination: this.inclination,
				longitude: this.longitude,
				anomaly: 0
			});
			const cost = raisedOrbit.velP - this.velP;
			return {
				"orbit": raisedOrbit,
				"valueOf": function(){return cost}
			}
		}
		else if(r < this.P){
			let raisedOrbit = new Orbit({
				system: this.system,
				periapsis: r,
				apoapsis: this.apoapsis,
				inclination: this.inclination,
				longitude: this.longitude,
				anomaly: 0
			});
			const cost = this.velA - raisedOrbit.velA;
			return {
				"orbit": raisedOrbit,
				"valueOf": function(){return cost}
			}
		}
	}
	circularize(r){
		let velocity = this.v(r);
		let circularOrbit = new Orbit({
			system: this.system,
			radius: r
		});
		return {
			"orbit": circularOrbit,
			"radial": velocity.radial,
			"tangential": circularOrbit.vela - velocity.tangential,
			"valueOf": function(){return Math.hypot(velocity.radial,circularOrbit.vela - velocity.tangential)}
		}
	}
	static synodic(first,second){
		let t1;
		let t2;
		if(first instanceof System){
			t1 = first.orbit.period
		}
		else if(first instanceof Orbit){
			t1 = first.period
		}
		else{
			t1 = first
		}
		if(second instanceof System){
			t2 = second.orbit.period
		}
		else if(second instanceof Orbit){
			t2 = second.period
		}
		else{
			t2 = second
		}
		return new Time(Math.abs(1/(1/t1 - 1/t2)))
	}
}

[
{
	name: "UNIT",
	type: "mathematical object",
	GM: 1,
	radius: 0
},
{
	name: "NULL",
	type: "mathematical object",
	GM: 0,
	radius: 0
},
{
	name: "Sun",
	type: "star",
	GM: 1.32712440018e20,
	mass: 1.98855e30,
	radiusEquator: 695700e3,
	area: 6.09e18,
	volume: 1.41e27,
	density: 1.408e3,
	gravity: 274,
	color: "#FFFF00",
	satellites: ["Mercury","Venus","Earth","Mars","Ceres","Vesta","Jupiter","Saturn","Uranus","Neptune","Pluto","Eris"]
},
{
	name: "Earth",
	type: "terrestrial planet",
	GM: 3.986004418e14,
	mass: 5.97237e24,
	radius: 6371.0e3,
	radiusPolar: 6356.8e3,
	radiusEquator: 6378.1e3,
	volume: 1.08321e21,
	density: 5.514e3,
	area: 510072000e6,
	gravity: 9.807,
	period: 0.99726968 * 86400,
	lowOrbitAltitude: 200e3,
	satellites: ["Moon"],
	orbit: {
		system: "Sun",
		apoapsis: 152100000e3,
		periapsis: 147095000e3,
		eccentricity: 0.0167086,
		period: 31558149.7635
	},
	geology: {
		liquidCore: true,
		k2: 0.29525
	},
	color: "#2f6a69",
	atmosphere: {
		height: 100e3,
		layers:[
			{
				"name":"troposphere",
				"base":-611,
				"temperature":19,
				"lapseRate":6.5,
				"pressure":108900
			},
			{
				"name":"tropopause",
				"base":11019,
				"temperature":-56.5,
				"lapseRate":0,
				"pressure":22632
			},
			{
				"name":"stratosphere",
				"base":20063,
				"temperature":-56.5,
				"lapseRate":-1,
				"pressure":5474.9
			},
			{
				"name":"stratosphere2",
				"base":32162,
				"temperature":-44.5,
				"lapseRate":-2.8,
				"pressure":868.02
			},
			{
				"name":"stratopause",
				"base":47350,
				"temperature":-2.5,
				"lapseRate":0,
				"pressure":110.91
			},
			{
				"name":"mesosphere",
				"base":51413,
				"temperature":-2.5,
				"lapseRate":2.8,
				"pressure":66.939
			},
			{
				"name":"mesosphere2",
				"base":71802,
				"temperature":-58.5,
				"lapseRate":2,
				"pressure":3.9564
			},
			{
				"name":"mesopause",
				"base":86000,
				"temperature":-86.28,
				"lapseRate":0,
				"pressure":0.3734
			}
		],
		composition: [
			["N2",0.78084],
			["O2",0.20946],
			["Ar",0.00934],
			["CO2",0.000407],
			["Ne",0.00001818],
			["He",0.00000524],
			["CH4",0.0000018],
			["Kr",0.00000114],
			["H2",0.00000055]
		].map(e => ({component: e[0],amount: e[1]})),
		molarMass: 0.02896
	}
},
{
	name: "Moon",
	type: "moon",
	GM: 4.9048695e12,
	mass: 7.342e22,
	radius: 1737.1e3,
	radiusPolar: 1736.0e3,
	radiusEquator: 1738.1e3,
	volume: 2.1958e19,
	density: 3.344e3,
	area: 3.793e13,
	gravity: 1.62,
	axialTilt: 6.687*(Math.PI/180),
	axialTiltEcliptic: 1.5424*(Math.PI/180),
	axialTiltParent: 24*(Math.PI/180),
	period: 27.321661*86400,
	tidalLock: true,
	geology: {
		k2: 0.02405
	},
	color: "#1F1F1F",
	orbit: {
		system: "Earth",
		apoapsis: 405400e3,
		periapsis: 362600e3,
		semiMajor: 384399e3,
		eccentricity: 0.0549,
		period: 27.321661*86400,
		inclination: 5.145*(Math.PI/180)
	}
},
{
	name: "Mercury",
	type: "terrestrial planet",
	GM: 2.2032e13,
	mass: 3.3011e23,
	radius: 2439.7e3,
	area: 7.48e13,
	volume: 6.083e19,
	density: 5.427e3,
	gravity: 3.7,
	axialTilt: 0.034*(Math.PI/180),
	period: 58.646 * 86400,
	color: "#1a1a1a",
	geology: {
		k2: 0.451
	},
	orbit: {
		system: "Sun",
		apoapsis: 69816900e3,
		periapsis: 46001200e3,
		semiMajor: 57909050e3,
		eccentricity: 0.205630,
		period: 87.9691*86400,
		inclination: 7.005*(Math.PI/180),
		longitude: 48.331*(Math.PI/180),
		argument: 29.124*(Math.PI/180),
	}
},
{
	name: "Venus",
	type: "terrestrial planet",
	GM: 3.24859e14,
	radius: 6051.8e3,
	period: -243.025 * 86400,
	geology: {
		k2: 0.295
	},
	color: "#e6e6e6",
	atmosphere: {
		layers: [
			{
				"name": "troposphere",
				"base": 0,
				"temperature": 462,
				"pressure": 92.1*101325
			},
			{
				"base": 5000,
				"temperature": 424,
				"pressure": 66.65*101325
			},
			{
				"base": 10000,
				"temperature": 385,
				"pressure": 47.39*101325
			},
			{
				"base": 15000,
				"temperature": 348,
				"pressure": 33.04*101325
			},
			{
				"base": 20000,
				"temperature": 306,
				"pressure": 22.52*101325
			},
			{
				"base": 25000,
				"temperature": 264,
				"pressure": 14.93*101325
			},
			{
				"base": 30000,
				"temperature": 222,
				"pressure": 9.851*101325
			},
			{
				"base": 35000,
				"temperature": 180,
				"pressure": 5.917*101325
			},
			{
				"base": 40000,
				"temperature": 143,
				"pressure": 3.501*101325
			},
			{
				"base": 45000,
				"temperature": 110,
				"pressure": 1.979*101325
			},
			{
				"base": 50000,
				"temperature": 75,
				"pressure": 1.066*101325
			},
			{
				"base": 55000,
				"temperature": 27,
				"pressure": 0.5314*101325
			},
			{
				"base": 60000,
				"temperature": -10,
				"pressure": 0.2357*101325
			},
			{
				"base": 65000,
				"temperature": -30,
				"lapseRate": 8.4,
				"pressure": 0.09765*101325
			},
			{
				"base": 70000,
				"temperature": -43,
				"pressure": 0.0369*101325
			},
			{
				"base": 80000,
				"temperature": -76,
				"pressure": 0.00476*101325
			},
			{
				"base": 90000,
				"temperature": -104,
				"pressure": 0.0003736*101325
			},
			{
				"base": 100000,
				"temperature": -112,
				"pressure": 0.0000266*101325
			}
		],
		height: 220000,
		composition: [
			["CO2",0.965],
			["N2",0.035],
			["SO2",0.000150],
			["Ar",0.000070],
			["H2O",0.000020],
			["CO",0.000017],
			["He",0.000012],
			["Ne",0.000007],
			["HCl",0.00000035],
			["HF",0.000000003]
		].map(e => ({component: e[0],amount: e[1]}))
	},
	orbit: {
		system: "Sun",
		apoapsis: 108939000e3,
		periapsis: 107477000e3,
		semiMajor: 108208000e3,
		eccentricity: 0.006772
	}
},
{
	name: "Mars",
	type: "terrestrial planet",
	GM: 4.282837e13,
	mass: 6.4171e23,
	radius: 3389.5e3,
	radiusPolar: 3376.2e3,
	radiusEquator: 3396.2e3,
	volume: 1.6318e20,
	density: 3.9335e3,
	area: 144798500e6,
	gravity: 3.72076,
	lowOrbitAltitude: 200e3,
	axialTilt: 25.19*(Math.PI/180),
	period: 1.025957*86400,
	satellites: ["Phobos","Deimos"],
	geology: {
		k2: 0.173
	},
	color: "#993d00",
	atmosphere: {
		height: 100000,
		layers: [
			{
				"name": "troposphere",
				"base": 0 ,
				"temperature": -58,
				"lapseRate": 2.5,
				"pressure": 610
			},
		],
		composition: [
			["CO2",0.949],
			["N2",0.026],
			["Ar",0.019],
			["O2",0.00174],
			["CO",0.000747],
			["H2O",0.0003]
		].map(e => ({component: e[0],amount: e[1]})),
	},
	orbit: {
		system: "Sun",
		apoapsis: 249200000e3,
		periapsis: 206700000e3,
		semiMajor: 227939200e3,
		eccentricity: 0.0934,
		period: 686.971*86400,
		inclination: 1.850*(Math.PI/180)
	}
},
{
	name: "Phobos",
	type: "moon",
	radius: 11.2667e3,
	axis: [27e3,22e3,18e3],
	orbit: {
		system: "Mars",
		periapsis: 9234.42e3,
		apoapsis: 9517.58e3,
		semiMajor: 9376e3,
		eccentricity: 0.0151
	}
},
{
	name: "Deimos",
	type: "moon",
	orbit: {
		system: "Mars",
		periapsis: 23455.5e3,
		apoapsis: 23470.9e3,
		semiMajor: 23463.2e3,
		eccentricity: 0.00033,
		period: 109123.2
	}
},
{
	name: "Ceres",
	type: "dwarf planet",
	period: 9.074170 * 3600,
	radius: 469.73e3,
	GM: 6.26325e10,
	mass: 9.3835e20,
	area: 2770000e6,
	density: 2.162e3,
	
	orbit: {
		system: "Sun",
		apoapsis: 445749000e3,
		periapsis: 382774000e3,
		semiMajor: 414261000e3,
		eccentricity: 0.07600902910,
		period: 466.6*86400,
		inclination: 10.59406704*(Math.PI/180),
		longitude: 80.3055316*(Math.PI/180),
		argument: 73.5976941*(Math.PI/180)
	}
},
{
	name: "Vesta",
	type: "asteroid",
	orbit: {
		system: "Sun"
	}
},
{
	name: "Jupiter",
	type: "gas planet",
	GM: 1.26686534e17,
	mass: 1.8982e27,
	radius: 69911e3,
	radiusEquator: 71492e3,
	radiusPolar: 66854e3,
	area: 6.1419e16,
	volume: 1.4313e24,
	density: 1.326e3,
	gravity: 24.79,
	period: 9.925 * 3600,
	axialTilt: 3.13*(Math.PI/180),
	satelittes: ["Io","Europa","Ganymede","Callisto"],
	color: "#b07f35",
	orbit: {
		system: "Sun",
		apoapsis: 816.62e9,
		periapsis: 740.52e9,
		semiMajor: 778.57e9,
		eccentricity: 0.0489,
		inclination: 1.303*(Math.PI/180),
		longitude: 100.464*(Math.PI/180),
		argument: 273.867*(Math.PI/180)
	}
},
{
	name: "Saturn",
	type: "gas planet",
	GM: 3.7931187e16,
	radius: 58232e3,
	radiusPolar: 54364e3,
	radiusEquator: 60268e3,
	satellites: ["Titan"],
	geology: {
		k2: 0.390
	},
	color: "#b08f36",
	orbit: {
		system: "Sun",
		semiMajor: 1433.53e9,
		apoapsis: 1514.50e9,
		periapsis: 1352.55e9
	}
},
{
	name: "Uranus",
	type: "gas planet",
	GM: 5.793939e15,
	color: "#5580aa",
	radius: 25362e3,
	radiusPolar: 24973e3,
	radiusequator: 25559e3,
	orbit: {
		system: "Sun",
		period: 30687.153*86400,
		semiMajor: 2875.04e9,
		periapsis: 2742e9,
		apoapsis: 3008e9,
		eccentricity: 0.046381,
		inclination: 0.773,
		argument: 96.998857
	}
},
{
	name: "Neptune",
	type: "gas planet",
	GM: 6.836529e15,
	radius: 24622e3,
	satellites: ["Triton"],
	color: "#366896",
	orbit: {
		system: "Sun",
		period: 60190.03*86400,
		periapsis: 4459500000000,
		apoapsis: 4537300000000,
		semiMajor: 4504400000000,
		eccentricity: 0.009456,
	}
},
{
	name: "Pluto",
	type: "dwarf planet",
	GM: 8.71e11,
	mass: 1.303e22,
	radius: 1188.3e3,
	density: 1.854e3,
	satellites: ["Charon","Nix","Hydra"],
	atmosphere: {
		layers: [
			{
				"base": 0 ,
				"temperature": -229,
				"pressure": 1
			},
		]
	},
	orbit: {
		system: "Sun",
		apoapsis: 7.37593e12,
		periapsis: 4.43682e12,
		semiMajor: 5.90638e12,
		eccentricity: 0.2488,
		period: 90560*86400,
		inclination: 17.16*(Math.PI/180),
		longitude: 110.299*(Math.PI/180),
		argument: 113.834*(Math.PI/180),
	}
},
{
	name: "Eris",
	type: "dwarf planet",
	GM: 1.108e12,
	orbit: {
		system: "Sun",
		apoapsis: 14.579e12,
		periapsis: 5.725e12
	}
},
{
	name: "Io",
	type: "moon",
	mass: 8.931938e22,
	radius: 1821.6e3,
	tidalLock: true,
	geology: {
		k2: 0.015
	},
	orbit: {
		system: "Jupiter",
		apoapsis: 423400e3,
		periapsis: 420000e3,
		semiMajor: 423400e3,
		eccentricity: 0.0041
	}
},
{
	name: "Europa",
	type: "moon",
	mass: 4.799844e22,
	radius: 1560.8e3,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Jupiter",
		semiMajor: 670900e3,
		eccentricity: 0.009
	}
},
{
	name: "Ganymede",
	type: "moon",
	mass: 1.4819e23,
	radius: 2634.1e3,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Jupiter",
		semiMajor: 1070400e3,
		eccentricity: 0.0013
	}
},
{
	name: "Callisto",
	type: "moon",
	radius: 2410.3e3,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Jupiter",
		semiMajor: 1882700e3,
		eccentricity: 0.0074
	}
},
{
	name: "Titan",
	type: "moon",
	radius: 2574.73e3,
	mass: 1.3452e23,
	tidalLock: true,
	geology: {},
	atmosphere: {
		layers: [],
		height: 1200000
	},
	orbit: {
		system: "Saturn",
		semiMajor: 1221870e3,
		eccentricity: 0.0288,
		inclination: 0.34854,
	}
},
{
	name: "Enceladus",
	type: "moon",
	radius: 252.1e3,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Saturn",
		semiMajor: 237948e3,
		eccentricity: 0.0047
	}
},
{
	name: "Mimas",
	type: "moon",
	radius: 198.2e3,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Saturn",
		semiMajor: 185539e3,
		eccentricity: 0.0196
	}
},
{
	name: "Tethys",
	type: "moon",
	radius: 531.1e3,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Saturn",
		semiMajor: 294619e3,
		eccentricity: 0.0001
	}
},
{
	name: "Dione",
	type: "moon",
	radius: 561.4e3,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Saturn",
		semiMajor: 377396e3,
		eccentricity: 0.0022
	}
},
{
	name: "Rhea",
	type: "moon",
	mass: 2.305518e21,
	radius: 763.8e3,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Saturn",
		semiMajor: 527108e3,
		eccentricity: 0.0012583
	}
},
{
	name: "Iapetus",
	type: "moon",
	mass: 1.805635e21,
	radius: 734.5e3,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Saturn",
		semiMajor: 3560820e3,
		eccentricity: 0.0276812
	}
},
{
	name: "Ariel",
	type: "moon",
	radius: 578.9e3,
	mass: 1.353e21,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Uranus",
		semiMajor: 191020e3,
		eccentricity: 0.0012
	}
},
{
	name: "Umbriel",
	type: "moon",
	radius: 584.7e3,
	mass: 1.172e21,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Uranus",
		semiMajor: 266000e3,
		eccentricity: 0.0039
	}
},
{
	name: "Titania",
	type: "moon",
	radius: 788.4e3,
	mass: 3.527e21,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Uranus",
		semiMajor: 435910e3,
		eccentricity: 0.0011
	}
},
{
	name: "Oberon",
	type: "moon",
	radius: 761.4e3,
	mass: 3.014e21,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Uranus",
		semiMajor: 583520e3,
		eccentricity: 0.0014
	}
},
{
	name: "Miranda",
	type: "moon",
	radius: 235.8e3,
	mass: 6.59e19,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Uranus",
		semiMajor: 129390e3,
		eccentricity: 0.0013
	}
},
{
	name: "Triton",
	type: "moon",
	radius: 1353.4e3,
	mass: 2.14e22,
	tidalLock: true,
	geology: {},
	orbit: {
		system: "Neptune",
		semiMajor: 354759e3,
		eccentricity: 0.000016
	}
}
].forEach(thing => {
	let system = new System(thing);
	if(thing.orbit){
		try{
			thing.orbit = new Orbit(thing.orbit)
		}
		catch(e){
			if(e === "invalid orbit"){
				console.log("orbit of",thing.name,"failed checks");
			}
			else{
				console.log("unknown error for orbit of",thing.name);
			}
			throw "orbit error"
		}
	}
	if(thing.atmosphere){
		thing.atmosphere = new Atmosphere(thing.atmosphere);
	}
	thing.geology = new Geology(thing.geology,system);
	systems.set(thing.name,system);
});

let earth = systems.get("Earth");
let moon = systems.get("Moon");
let sun = systems.get("Sun");
let venus = systems.get("Venus");
let mercury = systems.get("Mercury");
let mars = systems.get("Mars");
let ceres = systems.get("Ceres");
let jupiter = systems.get("Jupiter");
let saturn = systems.get("Saturn");
let uranus = systems.get("Uranus");
let neptune = systems.get("Neptune");

class Orbiter{
	constructor(properties){
		this.x = properties.x || 0;
		this.y = properties.y || 0;
		this.z = properties.z || 0;
		this.vx = properties.vx || 0;
		this.vy = properties.vy || 0;
		this.vz = properties.vz || 0;
	}
	static stepSizeCompleter(config){
		config.time = config.time || (config.stepSize || 1) * (config.steps || 1);
		config.steps = config.steps || Math.ceil(config.time / (config.stepSize || 1));
		config.stepSize = config.time / config.steps;
		return config
	}
	progress_1body(config){
		this.x += this.vx * config.time;
		this.y += this.vy * config.time;
		this.z += this.vz * config.time;
	}
	progress_2body(config){
	}
	progress_3body(config){
	}
	progress_nbody(config){
	}
	progress_patchedConics(config){
	}
	progress(config){
		config = Orbiter.stepSizeCompleter(config);
		if(config.mode === "1body"){
			this.progress_1body(config);
		}
		else if(config.mode === "2body"){
			this.progress_2body(config);
		}
		else if(config.mode === "3body"){
			this.progress_3body(config);
		}
		else if(config.mode === "nbody"){
			this.progress_nbody(config);
		}
		else if(config.mode === "patchedConics"){
			this.progress_patchedConics(config);
		}
	}
}

let yuri = function(angular,mu,surface,stationary){
	return (stationary*stationary - surface*surface)*angular*angular/2 + mu/stationary - mu/surface
}

let yuri2 = function(system){
	return yuri(2*Math.PI/system.period,system.GM,system.radius.equator,system.stationaryOrbit.a)
}

let integral = function(w,a,mu1,mu2,r){
	return w*w * (-r*r/2 + a*mu1*r/(mu1 + mu2)) + mu1/(r - a) - mu2/r
}

let yuriMoon = function(system,distance){
	let angular = 2*Math.PI/system.period;
	let surface = system.radius.equator;
	let semiMajor = system.orbit.a;
	let mu2 = system.GM;
	let mu1 = system.orbit.system.GM;
	return integral(angular,semiMajor,mu1,mu2,distance) - integral(angular,semiMajor,mu1,mu2,surface)
}

let asteroid_return = function(start_rad,peri){
	let marsPass = new Orbit({system: sun,apoapsis: start_rad, periapsis: peri}).inner(mars.orbit.a).velocity;
	let vinfinity = Math.hypot(marsPass.tangential - mars.orbit.vela, marsPass.radial);
	let turning = mars.turningAngle(mars.radius + 150000, vinfinity);
	let entryAngle = Math.acos(marsPass.radial/vinfinity);
	let exitAngle = entryAngle - turning;
	let v_t = Math.sin(exitAngle) * vinfinity + mars.orbit.vela;
	let v_r = - Math.cos(exitAngle) * vinfinity;
	let ret = Orbit.vectorsToOrbit(sun,[mars.orbit.a,0,0],[v_r,v_t,0]);

	console.log("burn1: ",- new Orbit({system: sun,apoapsis: start_rad, periapsis: peri}).velA + new Orbit({system: sun,apoapsis: start_rad, periapsis: start_rad}).velA);
	console.log("burn2: ",- new Orbit({system: sun,apoapsis: start_rad, periapsis: earth.orbit.a}).velA + new Orbit({system: sun,apoapsis: start_rad, periapsis: start_rad}).velA	);

	return ret;
}
/*

a(r) = w^2 * (a*mu1/(mu1 + mu2) - r) + mu2/r² - mu1/(a - r)^2

 w^2 * (a*mu1/(mu1 + mu2) - r) + mu2/r² = mu1/(a - r)^2

w² * (-x^2 / 2 + a*mu1*x/(mu1 + mu2)) + mu1/(x - a) - mu2/x

*/

