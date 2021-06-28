function tetherTool(DOMnode){

//part 1: utility functions

function create(
	HTMLtag,//<string>: The kind of DOM element to create
	classes,//(optional) <string>: The css class to give the element  OR   [<strings>]: A list of multiple class names
	//(the first string of the list could optionally start with a "#", in which case it will be an id instead)
	text,//(optional) <string>: The innerText to give the element
	appendLocation,//(optional) DOMnode: a node to immediately append the created element to
	cssText//(optional) <string>: Inline CSS to appy to the element
){
	let element = document.createElement(HTMLtag);
	if(Array.isArray(classes)){
		element.classList.add(...classes);
		if(classes.includes("newTab")){
			element.setAttribute("target","_blank")
		}
	}
	else if(classes){
		if(classes[0] === "#"){
			element.id = classes.substring(1)
		}
		else{
			element.classList.add(classes);
			if(classes === "newTab"){
				element.setAttribute("target","_blank")
			}
		}
	};
	if(text || text === 0){
		element.innerText = text
	};
	if(appendLocation && appendLocation.appendChild){
		appendLocation.appendChild(element)
	};
	if(cssText){
		element.style.cssText = cssText
	};
	return element
};

//part 2: css

let tetherToolStyle = create("style");
tetherToolStyle.id = "tetherToolCSS";
tetherToolStyle.type = "text/css";
document.head.appendChild(tetherToolStyle);
tetherToolStyle.textContent = `
.tetherToolMain .container{
	display: inline-block;
	-webkit-border-radius: 10px;
	-moz-border-radius: 10px;
	border-radius: 10px;
	border-style: solid;
	border-width: 1px;
	border-color: black;
	padding: 10px;
	margin: 5px;
}
.tetherToolMain{
	background-color: #F6F6F0;
	-webkit-border-radius: 10px;
	-moz-border-radius: 10px;
	border-radius: 10px;
	box-shadow: 5px 5px 4px #000000;
}
.tetherToolMain .canvas{
	border-style: solid;
	border-width: 1px;
	border-color: black;
	-moz-border-radius: 4px;
	-webkit-border-radius: 4px;
	border-radius: 4px;
}
.tetherToolMain .draggable{
	cursor: move;
}
.tetherToolMain .tooltip{
	cursor: help;
}
.tetherToolMain .svgText {
	pointer-events: none;
	-webkit-touch-callout: none;
	-webkit-user-select: none;
	-khtml-user-select: none;
	-moz-user-select: none;
	-ms-user-select: none;
	-o-user-select: none;
	user-select: none;
}
.tetherToolMain #materialOverride{
	display: none;
}
.tetherToolMain #overrideMaterial:checked ~ #materialOverride{
	display: block;
}
`

//part 3: HTML

let mainContainer = create("div","tetherToolMain",false,DOMnode);
create("h1",false,"Tethers",mainContainer);
create("span",false,"System",mainContainer);
let selector = create("select",false,false,mainContainer);
	selector.name = "system";
create("h3",false,"Vertical tether",mainContainer);
let formElement = create("form","#input",false,mainContainer);
	let formContainer = create("div","container",false,formElement);
		
	let warnings = create("p","#warnings",false,formElement,"color: red");


//part 4: planetary data

Array.from(systems.entries()).map(a => a[1]).filter(system => 
	system.properties.type !== "mathematical object"
	&& system.radius
	&& system.GM
).sort().forEach((system,index) => {
	var opt = create("option",false,system,selector);
	if(system.properties.name === "Moon"){//The Moon should be default
		selector.selectedIndex = index
	}
})

//part 5: tether tool

const svgNS = "http://www.w3.org/2000/svg";
var illustration = document.getElementById("svg_ilu");
var form = formElement.elements;
var moon = systems.get(selector.options[selector.selectedIndex].value);//parent object of the tether

var tether = {//remember to also update the select list in the HTML
	materials: {
		zylon:   [5800000000,1540],
		aramid:  [3620000000,1440],
		steel:   [2617000000,8000],
		hppe:	[2400000000,970],
		mwcnt:   [30000000000,1350],
		xxmwcnt: [100000000000,1350]
	}
};
function calc(){
	moon = systems.get(selector.options[selector.selectedIndex].value);//parent object of the tether
	var warnings = "";//message the user about malformed input
	var gnum = function(field){//form parser
		return Number(form[field].value);
	};
	var out = function(id,message){//return output to containers
		document.getElementById(id).innerText = message;
	};
	var myRound = function(num,digits){
		return Math.floor(num*Math.pow(10,digits))/Math.pow(10,digits);
	};
//construction of tether object
	tether.foot = gnum("foot")*1000;
	tether.top = gnum("top")*1000;
	tether.centre = gnum("centre")*1000;
	tether.Length = tether.top - tether.foot;
	if(tether.foot < 0){//some checks about malformed input
		warnings  += "The tether foot is inside the Body!<br>";
	};
	if(tether.foot > tether.top){
		warnings  += "The tether foot can not be higher than the top!<br>";
	}
	else if(tether.centre < tether.foot || tether.centre > tether.top){
		warnings += "The centre of the tether must be between the endpoints!<br>";
	};
	tether.time = 2*Math.PI*Math.sqrt(Math.pow(tether.centre + moon.radius,3)/moon.GM);
	tether.centreVel = Math.sqrt(moon.GM/(moon.radius + tether.centre));
	tether.footVel = (tether.foot + moon.radius) * tether.centreVel/(tether.centre + moon.radius);
	tether.topVel = (tether.top + moon.radius) * tether.centreVel/(tether.centre + moon.radius);
	tether.angularVelocity = tether.centreVel/(tether.centre + moon.radius);

	tether.topOrbit = Math.sqrt(moon.GM/(tether.top + moon.radius));
	tether.footOrbit = Math.sqrt(moon.GM/(tether.foot + moon.radius));
	tether.escapeFromTop = Math.sqrt(
		tether.topVel * tether.topVel
		- tether.topOrbit * tether.topOrbit * 2
	);

	tether.hitTheGround = Math.sqrt(
		tether.footVel * tether.footVel
		+ 2*moon.GM / moon.radius
		- tether.footOrbit * tether.footOrbit*2
	);

	tether.newApoapsis = Math.pow(
		(tether.top + moon.radius) * tether.topVel
		,2
	)/(
		2*moon.GM
		- (tether.top + moon.radius) * tether.topVel * tether.topVel
	);
	tether.newPeriapsis = Math.pow(
		(tether.foot + moon.radius) * tether.footVel
		,2
	)/(
		2*moon.GM
		- (tether.foot + moon.radius) * tether.footVel * tether.footVel
	);
	tether.collisionAngle = Math.acos(
		tether.angularVelocity*(tether.foot + moon.radius)*(tether.foot/moon.radius + 1)/tether.hitTheGround
	);

	if(form.overrideMaterial.checked){
		tether.strength = (gnum("tensileStrength")/gnum("safety"))/gnum("density");
	}
	else{
		tether.strength = (tether.materials[form["material"].value][0]/gnum("safety")) / tether.materials[form["material"].value][1];
	};

	tether.lowAcceleration = moon.GM/Math.pow(moon.radius + tether.foot,2) - tether.footVel*tether.footVel/(moon.radius + tether.foot);
	tether.highAcceleration = - moon.GM/Math.pow(moon.radius + tether.top,2) + tether.topVel*tether.topVel/(moon.radius + tether.top);

	//acceleration = mu/r² - r * v_c²/r_c²
	//integral = -mu/r - r²v_c²/(2r_c²)
	//upper part formula = -mu/r_c - v_c²/2 + mu/x + x²v_c²/(2r_c²)
	//upper taper formula = e^((-mu/r_c - v_c²/2 + mu/x + x²v_c²/(2r_c²))/strength)
	//Integral[e^((a+b/x + c*x^2)/d),x]
	// or Integral[e^(a+b/x + c*x^2),x]

	//alternatively:
	//acceleration = mu/r² - w²r
	//integral = -mu/r - w²r²/2

//calculation of taper ratio
	var integral = function(gm,angularVelocity,foot,top){
		return Math.abs(
			(
				-gm/top - angularVelocity*angularVelocity*top*top/2
			) - (
				-gm/foot - angularVelocity*angularVelocity*foot*foot/2
			)
		);
	};

	tether.lowIntegral = integral(
		moon.GM,
		tether.angularVelocity,
		tether.foot + moon.radius,
		tether.centre + moon.radius
	);
	tether.highIntegral = integral(
		moon.GM,
		tether.angularVelocity,
		tether.centre + moon.radius,
		tether.top + moon.radius
	);
	tether.lowRatio = Math.pow(Math.E,tether.lowIntegral/tether.strength);
	tether.highRatio = Math.pow(Math.E,tether.highIntegral/tether.strength);

	var lowIteratorSum = 0;
	var highIteratorSum = 0;
	var iteratorLimit = 1000;
	var lowCrosses = [];
	var highCrosses = [];
//numverical way to find the tether mass. Also used for plotting the width of it
	for(var i=0;i<iteratorLimit;i++){
		var lowCross = Math.pow(
			Math.E,
			integral(moon.GM,tether.angularVelocity,tether.foot + moon.radius,
				tether.foot + i*(tether.centre - tether.foot)/iteratorLimit + moon.radius
			)
			/tether.strength
		);
		lowIteratorSum += lowCross;
		var highCross = Math.pow(
			Math.E,
			integral(moon.GM,tether.angularVelocity,tether.centre + moon.radius,
				tether.centre + i*(tether.top - tether.centre)/iteratorLimit + moon.radius
			)
			/tether.strength
		);
		highIteratorSum += highCross;
		if(i % 10 == 0){
			lowCrosses.push(lowCross);
			highCrosses.push(highCross);
		};
	};

	tether.lowMass = tether.lowAcceleration*(tether.centre-tether.foot)*(lowIteratorSum/iteratorLimit)/tether.strength;
	tether.highMass = tether.highAcceleration*(tether.top-tether.centre)*(highIteratorSum/iteratorLimit)/tether.strength;

//parsing of planets.json to find bodies reachable after release from the tether
	var potentialTargets = {siblings: [], aunts: []};
	if(moon.orbit){
		potentialTargets.parent = moon.orbit.system;
		potentialTargets.siblings = moon.orbit.system.satellites.filter(sat => sat.properties.name !== moon.properties.name) || [];
		if(moon.orbit.system.orbit){
			potentialTargets.aunts = moon.orbit.system.orbit.system.satellites.filter(sat => sat.properties.name !== potentialTargets.parent.properties.name)
		}
	};
	var targets = [];
	if(potentialTargets.parent){
		targets.push(
			{
				name : potentialTargets.parent.properties.name,
				vinf : moon.orbit.vela - Math.sqrt(
					moon.orbit.system.GM * (2/(moon.orbit.a) - 2/(moon.orbit.a + moon.orbit.system.radius))
				)
			}
		);
	};
	for(var i=0;i<potentialTargets.siblings.length;i++){
		targets.push(
			{
				name : potentialTargets.siblings[i].properties.name,
				vinf : Math.abs(moon.orbit.vela - Math.sqrt(
					moon.orbit.system.GM * (2/(moon.orbit.a) - 2/(moon.orbit.a + potentialTargets.siblings[i].orbit.a))
				))
			}
		);
	};
	for(var i=0;i<potentialTargets.aunts.length;i++){
		targets.push(
			{
				name : potentialTargets.aunts[i].properties.name,
				vinf : Math.sqrt(
					Math.pow(
						moon.orbit.system.orbit.vela
						- Math.sqrt(
							moon.orbit.system.orbit.system.GM * (2/(moon.orbit.system.orbit.a) - 2/(moon.orbit.system.orbit.a + potentialTargets.aunts[i].orbit.a))
						),
						2
					) + moon.orbit.vela*moon.orbit.vela*2
				) - moon.orbit.vela
			}
		);
	};

	var findRadiusFromVinf = function(vinf,angular){//a binary search
		var vinf_squared = vinf*vinf;
		var angular_squared = angular*angular;
		var val = moon.radius + tether.centre;
		var step = (moon.radius + tether.centre)/2;
		while(angular_squared*val*val - 2*moon.GM/val < vinf_squared){
			val += step;
			step *= 2;
		};
		for(var i=0;i<50;i++){
			step = step/2;
			if(angular_squared*val*val - 2*moon.GM/val < vinf_squared){
				val += step;
			}
			else{
				val -= step;
			};
		};
		return val;
	};
	while(document.getElementById("releases").hasChildNodes()){
		document.getElementById("releases").removeChild(document.getElementById("releases").lastChild);
	};
	for(var i=0;i<targets.length;i++){
		if(Number.isNaN(targets[i].vinf)){
			continue;
		};
		targets[i].location = findRadiusFromVinf(targets[i].vinf,tether.angularVelocity);
		var item = document.createElement("p");
		item.innerText = targets[i].name;
		item.innerText += " " + myRound((targets[i].location - moon.radius)/1000,2) + "km";
		document.getElementById("releases").appendChild(item);
	};

/*Some math:

(newAltitude + mr) * x * pay + central * (newAltitude + mr + centre - foot) * (newAltitude + mr + centre - foot) * x/(newAltitude + mr)
	= tether.dockedMomentum

x = tether.dockedMomentum/(
		(newAltitude + mr) * pay + central * (newAltitude + mr + centre - foot)²/(newAltitude + mr)
	)

tether.dockedEnergy =
	(
		tether.dockedMomentum/(
			(newAltitude + mr) * pay + central * (newAltitude + mr + centre - foot)²/(newAltitude + mr)
		)
	)² * pay
	+ (
		(newAltitude + mr + centre - foot) * tether.dockedMomentum/(
			(newAltitude + mr) * pay + central * (newAltitude + mr + centre - foot)²/(newAltitude + mr)
		)/(newAltitude + mr)
	)² * central
	- pay * moon.gm/(newAltitude + mr)
	- central * moon.gm/(newAltitude + mr + centre - foot)


*/

	out("lowAcceleration","Acceleration at tether foot: " + myRound(tether.lowAcceleration,5) + " m/s² (" + myRound(100*tether.lowAcceleration/moon.surfaceGravity,2) + "% of surface gravity)");
	out("highAcceleration","Acceleration at tether top: " + myRound(tether.highAcceleration,5) + " m/s² (" + myRound(100*tether.highAcceleration/moon.surfaceGravity,2) + "% of surface gravity)");
	out("period","Period: " + myRound(tether.time/3600,3) + " hours");
	out("footVel","Tether velocity at foot " + myRound(tether.footVel,3) + " m/s");
	out("centreVel","Tether velocity at centre " + myRound(tether.centreVel,3) + " m/s");
	out("topVel","Tether velocity at top " + myRound(tether.topVel,3) + " m/s");
	if(tether.newApoapsis > 0){
		out("escape","Does not reach escape velocity");
		out("newApoapsis","Apoapsis after release from tether: " + myRound(tether.newApoapsis/1000,3) + " km");
	}
	else{
		out("escape","Vinf after release from tether top " + myRound(tether.escapeFromTop,3) + " m/s (840 m/s is required to reach Earth from the Moon)");
		out("newApoapsis","Apoapsis after release from tether: " + myRound(tether.newApoapsis/1000,3) + " km (negative means escape)");
	};
	if(tether.newPeriapsis < moon.radius){
		out("hitTheGround","Velocity at surface after release from tether foot " + myRound(tether.hitTheGround,3) + " m/s (collision angle " + myRound(180*(Math.PI/2 -tether.collisionAngle)/Math.PI,2) + "º from vertical)");
	}
	else{
		out("hitTheGround","");
	};
	out("lowRatio","foot-centre taper ratio: " + myRound(tether.lowRatio,3));
	out("highRatio","centre-top taper ratio: " + myRound(tether.highRatio,3));
	out("lowMass","Mass of lower tether: "
		+ myRound(tether.lowMass,4)
		+ " x payload mass (anything that's not the tether itself is payload)");
	out("highMass","Mass of higher tether: "
		+ myRound(tether.highMass,4)
		+ " x payload mass (anything that's not the tether itself is payload)");
	out("lowClimber","Lower tether climber energy usage: " + myRound(tether.lowIntegral,0) + " J/kg");
	out("highClimber","Upper tether climber energy usage: " + myRound(tether.highIntegral,0) + " J/kg");

	out("warnings",warnings);

	while(illustration.lastChild){//clean svg
		illustration.removeChild(illustration.lastChild);
	};
	minimum_x_coordinate = Math.max(-moon.radius/900,(moon.radius - (tether.top + tether.foot)/2)/1000);
	minimum_y_coordinate = Math.max(-moon.radius/900,-tether.Length/4000);
	minimum_x_distance = Math.min((moon.radius + tether.top)/850,tether.top/300) + Math.min(tether.Length/4000,moon.radius/1000);
	minimum_y_distance = Math.min(moon.radius/450,tether.Length/2000);
	illustration.setAttributeNS(null,"viewBox",minimum_x_coordinate + " " + minimum_y_coordinate + " " + minimum_x_distance + " " + minimum_y_distance);


//lots of svg drawing:
	var addText = function(container,content,x,y,size,color){
		var newMarker = document.createElementNS(svgNS,"text");
		newMarker.setAttributeNS(null,"x",x);
		newMarker.setAttributeNS(null,"y",y);
		newMarker.setAttributeNS(null,"font-size",size);
		newMarker.setAttributeNS(null,"class","svgText");
		newMarker.setAttributeNS(null,"fill",color);

		var textNode = document.createTextNode(content);
		newMarker.appendChild(textNode);
		container.appendChild(newMarker);
	};
	var addToolTip = function(container,content){
		var newToolTip = document.createElementNS(svgNS,"title");
		container.setAttributeNS(null,"class",container.getAttributeNS(null,"class") + " tooltip");
		var textNode = document.createTextNode(content);
		newToolTip.appendChild(textNode);
		container.appendChild(newToolTip);
	};
	var svg_planet = document.createElementNS(svgNS,"g");
	var svg_planet_main = document.createElementNS(svgNS,"circle");
	var svg_planet_atmosphere = document.createElementNS(svgNS,"circle");

	svg_planet_main.setAttributeNS(null,"cx",0);
	svg_planet_main.setAttributeNS(null,"cy",0);
	svg_planet_atmosphere.setAttributeNS(null,"cx",0);
	svg_planet_atmosphere.setAttributeNS(null,"cy",0);
	svg_planet_main.setAttributeNS(null,"fill",moon.color);
	svg_planet_atmosphere.setAttributeNS(null,"fill","#E0E0E0");
	svg_planet_main.setAttributeNS(null,"r",moon.radius/1000);
	svg_planet_atmosphere.setAttributeNS(null,"r",(moon.radius + (moon.atmosphere ? moon.atmosphere.height : 0))/1000);
	addToolTip(svg_planet_main,selector.options[selector.selectedIndex].value);

	svg_planet.appendChild(svg_planet_atmosphere);
	svg_planet.appendChild(svg_planet_main);

	var svg_scale = tether.Length/100000;

	svg_tether = document.createElementNS(svgNS,"g");
	var svg_tether_main = document.createElementNS(svgNS,"line");
	svg_tether_main.setAttributeNS(null,"x1",(moon.radius + tether.foot)/1000);
	svg_tether_main.setAttributeNS(null,"y1",0);
	svg_tether_main.setAttributeNS(null,"x2",(moon.radius + tether.top)/1000);
	svg_tether_main.setAttributeNS(null,"y2",0);
	svg_tether_main.setAttributeNS(null,"stroke","black");
	svg_tether_main.setAttributeNS(null,"stroke-linecap","butt");
	svg_tether_main.setAttributeNS(null,"stroke-width",svg_scale);
	svg_tether.appendChild(svg_tether_main);

	marker_container = document.createElementNS(svgNS,"g");
	text_marker_container = document.createElementNS(svgNS,"g");

	addText(
		text_marker_container,
		"Foot",
		(tether.foot + moon.radius)/1000,
		svg_scale*3,
		svg_scale*2.5,
		"black"
	);
	addText(
		text_marker_container,
		"Top",
		(tether.top + moon.radius)/1000,
		svg_scale*3,
		svg_scale*2.5,
		"black"
	);
	addText(
		text_marker_container,
		"Anchor",
		(tether.centre + moon.radius)/1000,
		svg_scale*3,
		svg_scale*2.5,
		"black"
	);

	var foot_marker = document.createElementNS(svgNS,"circle");
	foot_marker.setAttributeNS(null,"cx",(tether.foot + moon.radius)/1000);
	foot_marker.setAttributeNS(null,"cy",0);
	foot_marker.setAttributeNS(null,"fill","red");
	foot_marker.setAttributeNS(null,"class","draggable");
	foot_marker.setAttributeNS(null,"id","foot");
	foot_marker.setAttributeNS(null,"onmousedown","selectElement(evt)");
	foot_marker.setAttributeNS(null,"r",svg_scale*1.4);

	var anchor_marker = document.createElementNS(svgNS,"circle");
	anchor_marker.setAttributeNS(null,"cx",(tether.centre + moon.radius)/1000);
	anchor_marker.setAttributeNS(null,"cy",0);
	anchor_marker.setAttributeNS(null,"fill","#900000");
	anchor_marker.setAttributeNS(null,"class","draggable");
	anchor_marker.setAttributeNS(null,"id","anchor");
	anchor_marker.setAttributeNS(null,"onmousedown","selectElement(evt)");
	anchor_marker.setAttributeNS(null,"r",svg_scale*1.4);

	var top_marker = document.createElementNS(svgNS,"circle");
	top_marker.setAttributeNS(null,"cx",(tether.top + moon.radius)/1000);
	top_marker.setAttributeNS(null,"cy",0);
	top_marker.setAttributeNS(null,"fill","red");
	top_marker.setAttributeNS(null,"class","draggable");
	top_marker.setAttributeNS(null,"id","top");
	top_marker.setAttributeNS(null,"onmousedown","selectElement(evt)");
	top_marker.setAttributeNS(null,"r",svg_scale*1.4);

	var descent_trajectory = document.createElementNS(svgNS,"path");
	var semiMajor = (tether.newPeriapsis + moon.radius + tether.foot)/2000;
	var semiMinor = Math.sqrt(
		semiMajor*semiMajor
		- Math.pow((tether.foot + moon.radius - tether.newPeriapsis)/2000,2)
	);
	var descent_path = "M" + (moon.radius + tether.foot)/1000 + " 0a" + semiMajor + " " + semiMinor + " 0 0 0 " + 2*-semiMajor + " 0a" + semiMajor + " " + semiMinor + " 0 0 0 " + semiMajor*2 + " 0";
	descent_trajectory.setAttributeNS(null,"d",descent_path);
	descent_trajectory.setAttributeNS(null,"stroke","green");
	descent_trajectory.setAttributeNS(null,"fill","none");
	descent_trajectory.setAttributeNS(null,"stroke-width",svg_scale/3);
	descent_trajectory.setAttributeNS(null,"stroke-dasharray",svg_scale/2);
	if(tether.newPeriapsis < moon.radius){
		var collision_point = document.createElementNS(svgNS,"circle");
		var angle = Math.asin(
			(2*(tether.foot + moon.radius)*tether.newPeriapsis - moon.radius*(tether.foot + moon.radius + tether.newPeriapsis))
			/(moon.radius*(tether.foot + moon.radius - tether.newPeriapsis))
		) + Math.PI/2;
		var secondAngle = angle + tether.collisionAngle - Math.PI/2;
		collision_point.setAttributeNS(null,"cx",Math.cos(angle)*moon.radius/1000);
		collision_point.setAttributeNS(null,"cy",-Math.sin(angle)*moon.radius/1000);
		collision_point.setAttributeNS(null,"fill","blue");
		collision_point.setAttributeNS(null,"r",svg_scale);
		var collision_line = document.createElementNS(svgNS,"line");
		collision_line.setAttributeNS(null,"x1",0);
		collision_line.setAttributeNS(null,"y1",0);
		collision_line.setAttributeNS(null,"x2",Math.cos(angle)*moon.radius/1000);
		collision_line.setAttributeNS(null,"y2",-Math.sin(angle)*moon.radius/1000);
		collision_line.setAttributeNS(null,"fill","none");
		collision_line.setAttributeNS(null,"stroke","pink");
		collision_line.setAttributeNS(null,"stroke-linecap","round");
		collision_line.setAttributeNS(null,"stroke-width",svg_scale/2);
		var reference_line = document.createElementNS(svgNS,"line");
		reference_line.setAttributeNS(null,"x1",0);
		reference_line.setAttributeNS(null,"y1",0);
		reference_line.setAttributeNS(null,"x2",moon.radius/1000);
		reference_line.setAttributeNS(null,"y2",0);
		reference_line.setAttributeNS(null,"fill","none");
		reference_line.setAttributeNS(null,"stroke","pink");
		reference_line.setAttributeNS(null,"stroke-linecap","round");
		reference_line.setAttributeNS(null,"stroke-width",svg_scale/2);
		var reference_angle = document.createElementNS(svgNS,"path");
		reference_angle.setAttributeNS(null,"d","M "+moon.radius/3000+" 0A"+moon.radius/3000+" "+moon.radius/3000+" 0 0 0 "+Math.cos(angle)*moon.radius/3000+" "+Math.sin(angle)*moon.radius/-3000);
		reference_angle.setAttributeNS(null,"fill","none");
		reference_angle.setAttributeNS(null,"stroke","pink");
		reference_angle.setAttributeNS(null,"stroke-linecap","round");
		reference_angle.setAttributeNS(null,"stroke-width",svg_scale/2);
		var collision_angle = document.createElementNS(svgNS,"path");
		collision_angle.setAttributeNS(
			null,"d",
			"M "+2*Math.cos(angle)*moon.radius/3000+" "
			+2*-Math.sin(angle)*moon.radius/3000+"A"
			+moon.radius/3000+" "+moon.radius/3000+" 0 0 1 "
			+(Math.cos(angle)*moon.radius/1000 - moon.radius*Math.cos(secondAngle)/3000)+" "
			+(Math.sin(angle)*-moon.radius/1000 + moon.radius*Math.sin(secondAngle)/3000)
		);
		collision_angle.setAttributeNS(null,"fill","none");
		collision_angle.setAttributeNS(null,"stroke","red");
		collision_angle.setAttributeNS(null,"stroke-linecap","round");
		collision_angle.setAttributeNS(null,"stroke-width",svg_scale/2);
		var velocity_line = document.createElementNS(svgNS,"line");
		velocity_line.setAttributeNS(null,"x1",Math.cos(angle)*moon.radius/1000 - moon.radius*Math.cos(secondAngle)/1000);
		velocity_line.setAttributeNS(null,"y1",-Math.sin(angle)*moon.radius/1000 + moon.radius*Math.sin(secondAngle)/1000);
		velocity_line.setAttributeNS(null,"x2",Math.cos(angle)*moon.radius/1000 + moon.radius*Math.cos(secondAngle)/1000);
		velocity_line.setAttributeNS(null,"y2",-Math.sin(angle)*moon.radius/1000 - moon.radius*Math.sin(secondAngle)/1000);
		velocity_line.setAttributeNS(null,"fill","none");
		velocity_line.setAttributeNS(null,"stroke","red");
		velocity_line.setAttributeNS(null,"stroke-linecap","round");
		velocity_line.setAttributeNS(null,"stroke-width",svg_scale/2);
		addText(
			text_marker_container,
			myRound(180*angle/Math.PI,1) + "º",
			Math.cos(angle/2)*moon.radius/2500,
			-Math.sin(angle/2)*moon.radius/2500,
			svg_scale*2.5,
			"pink"
		);
		addText(
			text_marker_container,
			myRound(180*(Math.PI/2 -tether.collisionAngle)/Math.PI,1) + "º",
			Math.cos(angle)*moon.radius/1000 - moon.radius*Math.cos((secondAngle + angle)/2)/1800,
			-Math.sin(angle)*moon.radius/1000 + moon.radius*Math.sin((secondAngle + angle)/2)/1800,
			svg_scale*2.5,
			"red"
		);

		marker_container.appendChild(collision_angle);
		marker_container.appendChild(collision_line);
		marker_container.appendChild(reference_line);
		marker_container.appendChild(velocity_line);
		marker_container.appendChild(reference_angle);
		marker_container.appendChild(collision_point);
	};

	var anchor_trajectory = document.createElementNS(svgNS,"circle");
	anchor_trajectory.setAttributeNS(null,"stroke","green");
	anchor_trajectory.setAttributeNS(null,"fill","none");
	anchor_trajectory.setAttributeNS(null,"cx",0);
	anchor_trajectory.setAttributeNS(null,"cy",0);
	anchor_trajectory.setAttributeNS(null,"r",(tether.centre + moon.radius)/1000);
	anchor_trajectory.setAttributeNS(null,"stroke-width",svg_scale/3);
	anchor_trajectory.setAttributeNS(null,"stroke-dasharray",svg_scale/2);

	svg_tether.appendChild(descent_trajectory);
	svg_tether.appendChild(anchor_trajectory);
	marker_container.appendChild(foot_marker);
	marker_container.appendChild(anchor_marker);
	marker_container.appendChild(top_marker);

	for(var i=0;i<targets.length;i++){
		var marker = document.createElementNS(svgNS,"circle");
		marker.setAttributeNS(null,"cx",targets[i].location/1000);
		marker.setAttributeNS(null,"cy",0);
		marker.setAttributeNS(null,"fill","blue");
		marker.setAttributeNS(null,"r",svg_scale);
		addToolTip(marker,"Release point for " + targets[i].name + " transfer\nAltitude " + myRound((targets[i].location - moon.radius)/1000,2) + "km\nVinf " + myRound(targets[i].vinf,2) + "m/s");

		marker_container.appendChild(marker);

		addText(
			text_marker_container,
			targets[i].name + "↑",
			targets[i].location/1000,
			svg_scale*3,
			svg_scale*2.5,
			"blue"
		);
	};
	//construct path for tether width diagram
	var tether_width_scale = Math.min(
		Math.sqrt(tether.lowRatio)*svg_scale*4,
		Math.sqrt(tether.highRatio)*svg_scale*4,
		svg_scale*13
	);
	var lowerPath = "M" + (tether.foot + moon.radius)/1000 + " " + svg_scale*4.5 + "l" + (tether.centre - tether.foot)/1000 + " 0l0 " + tether_width_scale;
	var upperPath = "M" + (tether.top  + moon.radius)/1000 + " " + svg_scale*4.5 + "l" + (tether.centre - tether.top )/1000 + " 0l0 " + tether_width_scale;
	for(var i=1;i<lowCrosses.length;i++){
		lowerPath += "l-" + (tether.centre - tether.foot)/(1000*lowCrosses.length) + " " + (Math.sqrt(1/lowCrosses[i]) -  Math.sqrt(1/lowCrosses[i-1]))*tether_width_scale;
		upperPath += "l" + (tether.top - tether.centre)/(1000*highCrosses.length) + " " + (Math.sqrt(1/highCrosses[i]) -  Math.sqrt(1/highCrosses[i-1]))*tether_width_scale;
	};
	lowerPath += "l-" + (tether.centre - tether.foot)/(1000*lowCrosses.length) + " 0z";
	upperPath += "l" + (tether.top - tether.centre)/(1000*highCrosses.length) + " 0z";
	var tether_cross_low = document.createElementNS(svgNS,"path");
	tether_cross_low.setAttributeNS(null,"fill","black");
	tether_cross_low.setAttributeNS(null,"d",lowerPath);
	tether_cross_low.setAttributeNS(null,"id","lowerPath");
	tether_cross_low.setAttributeNS(null,"stroke","none");
	addToolTip(tether_cross_low,"Lower tether diameter\nDiameter ratio " + myRound(Math.sqrt(tether.lowRatio),3) + "\nCross section ratio " + myRound(tether.lowRatio,3));

	text_marker_container.appendChild(tether_cross_low);

	var tether_cross_high = document.createElementNS(svgNS,"path");
	tether_cross_high.setAttributeNS(null,"fill","black");
	tether_cross_high.setAttributeNS(null,"d",upperPath);
	tether_cross_high.setAttributeNS(null,"id","upperPath");
	tether_cross_high.setAttributeNS(null,"stroke","none");
	addToolTip(tether_cross_high,"Upper tether diameter\nDiameter ratio " + myRound(Math.sqrt(tether.highRatio),3) + "\nCross section ratio " + myRound(tether.highRatio,3));

	text_marker_container.appendChild(tether_cross_high);
	illustration.appendChild(svg_tether);
	illustration.appendChild(svg_planet);
	illustration.appendChild(marker_container);
	illustration.appendChild(text_marker_container);

};
var hardReload = function(){
	var animation = document.createElementNS(svgNS,"animate");
	animation.setAttributeNS(null,"attributeName","opacity");
	animation.setAttributeNS(null,"begin","indefinite");
	animation.setAttributeNS(null,"values","0;0;1");
	animation.setAttributeNS(null,"dur", 1.5);
	animation.setAttributeNS(null,"fill","freeze");

	var marker_animation = document.createElementNS(svgNS,"animate");
	marker_animation.setAttributeNS(null,"attributeName","opacity");
	marker_animation.setAttributeNS(null,"begin","indefinite");
	marker_animation.setAttributeNS(null,"values","0;0;0;0;1");
	marker_animation.setAttributeNS(null,"dur", 2);
	marker_animation.setAttributeNS(null,"fill","freeze");

	var text_marker_animation = document.createElementNS(svgNS,"animate");
	text_marker_animation.setAttributeNS(null,"attributeName","opacity");
	text_marker_animation.setAttributeNS(null,"begin","indefinite");
	text_marker_animation.setAttributeNS(null,"values","0;0;0;0;0;1");
	text_marker_animation.setAttributeNS(null,"dur", 2.5);
	text_marker_animation.setAttributeNS(null,"fill","freeze");

	var svg_rotate = document.createElementNS(svgNS,"animate");
	svg_rotate.setAttributeNS(null,"attributeName","viewBox");
	svg_rotate.setAttributeNS(null,"begin","indefinite");
	svg_rotate.setAttributeNS(null,"from",-moon.radius/900 + " " + (-moon.radius/900) + " " + (moon.radius*2 + tether.top)/225 + " " + moon.radius/112);
	svg_rotate.setAttributeNS(null,"to",minimum_x_coordinate + " " + minimum_y_coordinate + " " + minimum_x_distance + " " + minimum_y_distance);
	svg_rotate.setAttributeNS(null,"dur",0.5);

	svg_tether.appendChild(animation);
	text_marker_container.appendChild(text_marker_animation);
	marker_container.appendChild(marker_animation);
	illustration.appendChild(svg_rotate);

	svg_rotate.beginElement();
	marker_animation.beginElement();
	text_marker_animation.beginElement();
};
var selectedElement = 0;
function selectElement(evt){
	selectedElement = evt.target;
	illustration.setAttributeNS(null, "onmousemove", "moveElement(evt)");
	illustration.setAttributeNS(null, "onmouseup", "deselectElement(evt)");
};
function moveElement(evt){
	var pt = illustration.createSVGPoint(), svgP, circle;
	pt.x = evt.clientX;
	pt.y = evt.clientY;
	svgP = pt.matrixTransform(illustration.getScreenCTM().inverse());
	if(selectedElement.id == "anchor"){
		if(svgP.x > (tether.foot + moon.radius)/1000 && svgP.x < (tether.top + moon.radius)/1000){
			selectedElement.setAttributeNS(null,"cx",svgP.x);
			form["centre"].value = svgP.x - moon.radius/1000;
		};
	}
	else if(selectedElement.id == "foot"){
		if(svgP.x < (tether.centre + moon.radius)/1000 && svgP.x > moon.radius/1000){
			selectedElement.setAttributeNS(null,"cx",svgP.x);
			form["foot"].value = svgP.x - moon.radius/1000;
		};
	}
	else if(selectedElement.id == "top"){
		if(svgP.x > (tether.centre + moon.radius)/1000){
			selectedElement.setAttributeNS(null,"cx",svgP.x);
			form["top"].value = svgP.x - moon.radius/1000;
		};
	}
	else{//default behaviour of move
		selectedElement.setAttributeNS(null,"cx",svgP.x);
		selectedElement.setAttributeNS(null,"cy",svgP.y);
	};
};
function deselectElement(evt){
	if(selectedElement != 0){
		illustration.removeAttributeNS(null, "onmousemove");
		illustration.removeAttributeNS(null, "onmouseup");
		selectedElement = 0;
	};
	calc();
};
calc();hardReload();
selector.oninput = function(){calc();hardReload();}
}
