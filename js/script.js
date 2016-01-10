    
var scene;
var camera;
var renderer;
var projector;
var controls;
var light;
var onRenderFcts = [];
var Vector3 = THREE.Vector3
var planet
var showAirCurrents = false
var showPlateBoundaries = false
var showPlateMovements = false

var subdivisions = 40;
var distortionRate = 1;
var plateCount = 50;
var oceanicRate = 0.7;
var heatLevel = 0.9;
var moistureLevel = 1;

var extrudeLevel = 45;

var tileborder = 0.05

var rotate = false;
var rotationSpeed = 0.001


//
window.onload = function(){
    initializeScene();
    planet = generatePlanet(subdivisions,distortionRate,plateCount,oceanicRate,heatLevel,moistureLevel, extrudeLevel);
    drawPlanet(planet)
    renderScene(); 
    console.log(planet)     
}

function initializeScene(){
    // init renderer
    renderer    = new THREE.WebGLRenderer({
        antialias   : true
    });
    renderer.setClearColor(new THREE.Color('lightgrey'), 1)
    renderer.setFaceCulling(THREE.CullFaceFront, THREE.FrontFaceDirectionCW);
    renderer.sortObjects = false;

    renderer.setSize( window.innerWidth, window.innerHeight );
    document.body.appendChild( renderer.domElement );

    projector = new THREE.Projector();



    // init scene and camera
    scene = new THREE.Scene();
    
    camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.01,50000);
    camera.position.z =3000;

    controls = new THREE.OrbitControls(camera);
    addLight();
}

function addLight(){
    var ambientLight = new THREE.AmbientLight(0xFFFFFF);
    scene.add(ambientLight);
    
    directionalLight = new THREE.DirectionalLight(0xFFFFFF);
    directionalLight.position.set(-3, 3, 7).normalize();
    scene.add(directionalLight);
}

function renderScene(){
    
    // handle window resize
    window.addEventListener('resize', function(){
        renderer.setSize( window.innerWidth, window.innerHeight )
        camera.aspect   = window.innerWidth / window.innerHeight
        camera.updateProjectionMatrix()     
    }, false)
    // render the scene
    onRenderFcts.push(function(){
        renderer.render( scene, camera );       
    })

    // run the rendering loop
    var xAxis = new THREE.Vector3(1,1,1);
    var lastTimeMsec= null
    requestAnimationFrame(function animate(nowMsec){
        // keep looping
        requestAnimationFrame( animate );
        // measure time
        lastTimeMsec    = lastTimeMsec || nowMsec-1000/60
        var deltaMsec   = Math.min(200, nowMsec - lastTimeMsec)
        lastTimeMsec    = nowMsec
        // call each update function
        onRenderFcts.forEach(function(onRenderFct){
            onRenderFct(deltaMsec/1000, nowMsec/1000)
        })
        if(rotate)
        rotateAroundWorldAxis(planet.renderData.surface.renderObject, xAxis, rotationSpeed);;      
        
    })
}

// Rotate an object around an arbitrary axis in object space
var rotObjectMatrix;
function rotateAroundObjectAxis(object, axis, radians) {
    rotObjectMatrix = new THREE.Matrix4();
    rotObjectMatrix.makeRotationAxis(axis.normalize(), radians);
    object.matrix.multiply(rotObjectMatrix);
    object.rotation.setFromRotationMatrix(object.matrix);
}

var rotWorldMatrix;
// Rotate an object around an arbitrary axis in world space       
function rotateAroundWorldAxis(object, axis, radians) {
    rotWorldMatrix = new THREE.Matrix4();
    rotWorldMatrix.makeRotationAxis(axis.normalize(), radians);
    rotWorldMatrix.multiply(object.matrix);                // pre-multiply
    object.matrix = rotWorldMatrix;
    object.rotation.setFromRotationMatrix(object.matrix);
}


//PLANET
function generatePlanet(subDivisions,topologyDistortionRate, plateCount, oceanicRate,heatLevel,moistureLevel, extrudeLevel){
    console.log("Generating planet")
    var planet = {};
    var mesh;
    planet.material = [new THREE.MeshLambertMaterial({ color: new THREE.Color(0x000000), ambient: new THREE.Color(0xFFFFFF), vertexColors: THREE.VertexColors, })];

    //generate subdivided icosahedron
    console.log("Generating Mesh")
    planet.topology = generateSubdividedIcosahedron(subDivisions);

    console.log("Distorting Mesh")
    // distort mesh
    var totalDistortion = Math.ceil(planet.topology.edges.length * topologyDistortionRate);
    var remainingIterations = 6;
    while(remainingIterations>0){
        var iterationDistortion = Math.floor(totalDistortion / remainingIterations);
        totalDistortion -= iterationDistortion;
        distortMesh(planet.topology, iterationDistortion);
        relaxMesh(planet.topology, 0.5);
        --remainingIterations; 
    }

    console.log("Relaxing Mesh")
    var intervalIteration = 0
    var initialIntervalIteration = intervalIteration;
    
    var averageNodeRadius = Math.sqrt(4 * Math.PI / planet.topology.nodes.length);
    var minShiftDelta = averageNodeRadius / 50000 * planet.topology.nodes.length;
    var maxShiftDelta = averageNodeRadius / 50 * planet.topology.nodes.length;
    var priorShift;
    var currentShift = relaxMesh(planet.topology, 0.5);
    var x = 4;
    while(x>0){
        priorShift = currentShift;
        currentShift = relaxMesh(planet.topology, 0.5);
        var shiftDelta = Math.abs(currentShift - priorShift);
        if (shiftDelta >= minShiftDelta && intervalIteration - initialIntervalIteration < 300)
        {   
            x+=(maxShiftDelta - shiftDelta) / (maxShiftDelta - minShiftDelta)
        }
        x--  
    }
    
    // calculate centroids 
    console.log("Calculating Centroids")
    centroids(planet.topology)

    //reorder triangle nodes
    console.log("Reordering triangle nodes")
    reorderTriangleNodes(planet.topology)


    //compute dual polyhedron
    console.log("Generating topology")
    planet.tiledTopology = generatePlanetTiledTopology(planet.topology);

    
    console.log("Generating Planet Terrain")
    generatePlanetTerrain(planet, plateCount, oceanicRate);

    
    console.log("Generating Planet Weather")
    generatePlanetWeather(planet.tiledTopology,heatLevel,moistureLevel)
    console.log("Generating Planet Biomes")
    generatePlanetBiomes(planet.tiledTopology.tiles,1000)
    extrude(planet, extrudeLevel)

    planet.renderData = generatePlanetRenderData(planet.tiledTopology)


    return planet;
}

function drawPlanet(planet){
    setSurfaceRenderMode(planet,"terrain")
    var obj = planet.renderData.surface.renderObject
    obj.name = "planet"
    scene.add(obj)
}



//TOPOLOGY
function generateIcosahedron(){
    var phi = (1.0 + Math.sqrt(5.0)) / 2.0;
    var du = 1.0 / Math.sqrt(phi * phi + 1.0);
    var dv = phi * du;
    
    nodes =
    [
        { p: new Vector3(0, +dv, +du), e: [], f: [] },
        { p: new Vector3(0, +dv, -du), e: [], f: [] },
        { p: new Vector3(0, -dv, +du), e: [], f: [] },
        { p: new Vector3(0, -dv, -du), e: [], f: [] },
        { p: new Vector3(+du, 0, +dv), e: [], f: [] },
        { p: new Vector3(-du, 0, +dv), e: [], f: [] },
        { p: new Vector3(+du, 0, -dv), e: [], f: [] },
        { p: new Vector3(-du, 0, -dv), e: [], f: [] },
        { p: new Vector3(+dv, +du, 0), e: [], f: [] },
        { p: new Vector3(+dv, -du, 0), e: [], f: [] },
        { p: new Vector3(-dv, +du, 0), e: [], f: [] },
        { p: new Vector3(-dv, -du, 0), e: [], f: [] },
    ];
    
    edges =
    [
        { n: [  0,  1, ], f: [], },
        { n: [  0,  4, ], f: [], },
        { n: [  0,  5, ], f: [], },
        { n: [  0,  8, ], f: [], },
        { n: [  0, 10, ], f: [], },
        { n: [  1,  6, ], f: [], },
        { n: [  1,  7, ], f: [], },
        { n: [  1,  8, ], f: [], },
        { n: [  1, 10, ], f: [], },
        { n: [  2,  3, ], f: [], },
        { n: [  2,  4, ], f: [], },
        { n: [  2,  5, ], f: [], },
        { n: [  2,  9, ], f: [], },
        { n: [  2, 11, ], f: [], },
        { n: [  3,  6, ], f: [], },
        { n: [  3,  7, ], f: [], },
        { n: [  3,  9, ], f: [], },
        { n: [  3, 11, ], f: [], },
        { n: [  4,  5, ], f: [], },
        { n: [  4,  8, ], f: [], },
        { n: [  4,  9, ], f: [], },
        { n: [  5, 10, ], f: [], },
        { n: [  5, 11, ], f: [], },
        { n: [  6,  7, ], f: [], },
        { n: [  6,  8, ], f: [], },
        { n: [  6,  9, ], f: [], },
        { n: [  7, 10, ], f: [], },
        { n: [  7, 11, ], f: [], },
        { n: [  8,  9, ], f: [], },
        { n: [ 10, 11, ], f: [], },
    ];
    
    faces =
    [
        { n: [  0,  1,  8 ], e: [  0,  7,  3 ], },
        { n: [  0,  4,  5 ], e: [  1, 18,  2 ], },
        { n: [  0,  5, 10 ], e: [  2, 21,  4 ], },
        { n: [  0,  8,  4 ], e: [  3, 19,  1 ], },
        { n: [  0, 10,  1 ], e: [  4,  8,  0 ], },
        { n: [  1,  6,  8 ], e: [  5, 24,  7 ], },
        { n: [  1,  7,  6 ], e: [  6, 23,  5 ], },
        { n: [  1, 10,  7 ], e: [  8, 26,  6 ], },
        { n: [  2,  3, 11 ], e: [  9, 17, 13 ], },
        { n: [  2,  4,  9 ], e: [ 10, 20, 12 ], },
        { n: [  2,  5,  4 ], e: [ 11, 18, 10 ], },
        { n: [  2,  9,  3 ], e: [ 12, 16,  9 ], },
        { n: [  2, 11,  5 ], e: [ 13, 22, 11 ], },
        { n: [  3,  6,  7 ], e: [ 14, 23, 15 ], },
        { n: [  3,  7, 11 ], e: [ 15, 27, 17 ], },
        { n: [  3,  9,  6 ], e: [ 16, 25, 14 ], },
        { n: [  4,  8,  9 ], e: [ 19, 28, 20 ], },
        { n: [  5, 11, 10 ], e: [ 22, 29, 21 ], },
        { n: [  6,  9,  8 ], e: [ 25, 28, 24 ], },
        { n: [  7, 10, 11 ], e: [ 26, 29, 27 ], },
    ];
    
    for (var i = 0; i < edges.length; ++i)
        for (var j = 0; j < edges[i].n.length; ++j)
            nodes[j].e.push(i);
    
    for (var i = 0; i < faces.length; ++i)
        for (var j = 0; j < faces[i].n.length; ++j)
            nodes[j].f.push(i);
    
    for (var i = 0; i < faces.length; ++i)
        for (var j = 0; j < faces[i].e.length; ++j)
            edges[j].f.push(i);

    
    return { nodes: nodes, edges: edges, faces: faces };
}

function generateSubdividedIcosahedron(degree){
    var icosahedron = generateIcosahedron();
    
    var nodes = [];
    for (var i = 0; i < icosahedron.nodes.length; ++i){
        nodes.push({ p: icosahedron.nodes[i].p, e: [], f: [] });
    }
    
    var edges = [];
    for (var i = 0; i < icosahedron.edges.length; ++i){
        var edge = icosahedron.edges[i];
        edge.subdivided_n = [];
        edge.subdivided_e = [];
        var n0 = icosahedron.nodes[edge.n[0]];
        var n1 = icosahedron.nodes[edge.n[1]];
        var p0 = n0.p;
        var p1 = n1.p;
        var delta = p1.clone().sub(p0);
        nodes[edge.n[0]].e.push(edges.length);
        var priorNodeIndex = edge.n[0];
        for (var s = 1; s < degree; ++s)
        {
            var edgeIndex = edges.length;
            var nodeIndex = nodes.length;
            edge.subdivided_e.push(edgeIndex);
            edge.subdivided_n.push(nodeIndex);
            edges.push({ n: [ priorNodeIndex, nodeIndex ], f: [] });
            priorNodeIndex = nodeIndex;
            nodes.push({ p: slerp(p0, p1, s / degree), e: [ edgeIndex, edgeIndex + 1 ], f: [] });
        }
        edge.subdivided_e.push(edges.length);
        nodes[edge.n[1]].e.push(edges.length);
        edges.push({ n: [ priorNodeIndex, edge.n[1] ], f: [] });
    }

    var faces = [];
    for (var i = 0; i < icosahedron.faces.length; ++i){
        var face = icosahedron.faces[i];
        var edge0 = icosahedron.edges[face.e[0]];
        var edge1 = icosahedron.edges[face.e[1]];
        var edge2 = icosahedron.edges[face.e[2]];
        var point0 = icosahedron.nodes[face.n[0]].p;
        var point1 = icosahedron.nodes[face.n[1]].p;
        var point2 = icosahedron.nodes[face.n[2]].p;
        var delta = point1.clone().sub(point0);
        
        var getEdgeNode0 = (face.n[0] === edge0.n[0])
            ? function(k) { return edge0.subdivided_n[k]; }
            : function(k) { return edge0.subdivided_n[degree - 2 - k]; };
        var getEdgeNode1 = (face.n[1] === edge1.n[0])
            ? function(k) { return edge1.subdivided_n[k]; }
            : function(k) { return edge1.subdivided_n[degree - 2 - k]; };
        var getEdgeNode2 = (face.n[0] === edge2.n[0])
            ? function(k) { return edge2.subdivided_n[k]; }
            : function(k) { return edge2.subdivided_n[degree - 2 - k]; };

        var faceNodes = [];
        faceNodes.push(face.n[0]);
        for (var j = 0; j < edge0.subdivided_n.length; ++j)
            faceNodes.push(getEdgeNode0(j));
        faceNodes.push(face.n[1]);
        for (var s = 1; s < degree; ++s)
        {
            faceNodes.push(getEdgeNode2(s - 1));
            var p0 = nodes[getEdgeNode2(s - 1)].p;
            var p1 = nodes[getEdgeNode1(s - 1)].p;
            for (var t = 1; t < degree - s; ++t)
            {
                faceNodes.push(nodes.length);
                nodes.push({ p: slerp(p0, p1, t / (degree - s)), e: [], f: [], });
            }
            faceNodes.push(getEdgeNode1(s - 1));
        }
        faceNodes.push(face.n[2]);
        
        var getEdgeEdge0 = (face.n[0] === edge0.n[0])
            ? function(k) { return edge0.subdivided_e[k]; }
            : function(k) { return edge0.subdivided_e[degree - 1 - k]; };
        var getEdgeEdge1 = (face.n[1] === edge1.n[0])
            ? function(k) { return edge1.subdivided_e[k]; }
            : function(k) { return edge1.subdivided_e[degree - 1 - k]; };
        var getEdgeEdge2 = (face.n[0] === edge2.n[0])
            ? function(k) { return edge2.subdivided_e[k]; }
            : function(k) { return edge2.subdivided_e[degree - 1 - k]; };

        var faceEdges0 = [];
        for (var j = 0; j < degree; ++j)
            faceEdges0.push(getEdgeEdge0(j));
        var nodeIndex = degree + 1;
        for (var s = 1; s < degree; ++s)
        {
            for (var t = 0; t < degree - s; ++t)
            {
                faceEdges0.push(edges.length);
                var edge = { n: [ faceNodes[nodeIndex], faceNodes[nodeIndex + 1], ], f: [], };
                nodes[edge.n[0]].e.push(edges.length);
                nodes[edge.n[1]].e.push(edges.length);
                edges.push(edge);
                ++nodeIndex;
            }
            ++nodeIndex;
        }

        var faceEdges1 = [];
        nodeIndex = 1;
        for (var s = 0; s < degree; ++s)
        {
            for (var t = 1; t < degree - s; ++t)
            {
                faceEdges1.push(edges.length);
                var edge = { n: [ faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s], ], f: [], };
                nodes[edge.n[0]].e.push(edges.length);
                nodes[edge.n[1]].e.push(edges.length);
                edges.push(edge);
                ++nodeIndex;
            }
            faceEdges1.push(getEdgeEdge1(s));
            nodeIndex += 2;
        }

        var faceEdges2 = [];
        nodeIndex = 1;
        for (var s = 0; s < degree; ++s)
        {
            faceEdges2.push(getEdgeEdge2(s));
            for (var t = 1; t < degree - s; ++t)
            {
                faceEdges2.push(edges.length);
                var edge = { n: [ faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s + 1], ], f: [], };
                nodes[edge.n[0]].e.push(edges.length);
                nodes[edge.n[1]].e.push(edges.length);
                edges.push(edge);
                ++nodeIndex;
            }
            nodeIndex += 2;
        }
        
        nodeIndex = 0;
        edgeIndex = 0;
        for (var s = 0; s < degree; ++s)
        {
            for (t = 1; t < degree - s + 1; ++t)
            {
                var subFace = {
                    n: [ faceNodes[nodeIndex], faceNodes[nodeIndex + 1], faceNodes[nodeIndex + degree - s + 1], ],
                    e: [ faceEdges0[edgeIndex], faceEdges1[edgeIndex], faceEdges2[edgeIndex], ], };
                nodes[subFace.n[0]].f.push(faces.length);
                nodes[subFace.n[1]].f.push(faces.length);
                nodes[subFace.n[2]].f.push(faces.length);
                edges[subFace.e[0]].f.push(faces.length);
                edges[subFace.e[1]].f.push(faces.length);
                edges[subFace.e[2]].f.push(faces.length);
                faces.push(subFace);
                ++nodeIndex;
                ++edgeIndex;
            }
            ++nodeIndex;
        }
        
        nodeIndex = 1;
        edgeIndex = 0;
        for (var s = 1; s < degree; ++s)
        {
            for (t = 1; t < degree - s + 1; ++t)
            {
                var subFace = {
                    n: [ faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s + 2], faceNodes[nodeIndex + degree - s + 1], ],
                    e: [ faceEdges2[edgeIndex + 1], faceEdges0[edgeIndex + degree - s + 1], faceEdges1[edgeIndex], ], };
                nodes[subFace.n[0]].f.push(faces.length);
                nodes[subFace.n[1]].f.push(faces.length);
                nodes[subFace.n[2]].f.push(faces.length);
                edges[subFace.e[0]].f.push(faces.length);
                edges[subFace.e[1]].f.push(faces.length);
                edges[subFace.e[2]].f.push(faces.length);
                faces.push(subFace);
                ++nodeIndex;
                ++edgeIndex;
            }
            nodeIndex += 2;
            edgeIndex += 1;
        }
    }

    return { nodes: nodes, edges: edges, faces: faces };
}

function centroids(topology){
    //calculate face centroids
    for (var i = 0; i < topology.faces.length; ++i)
    {
        var face = topology.faces[i];
        var p0 = topology.nodes[face.n[0]].p;
        var p1 = topology.nodes[face.n[1]].p;
        var p2 = topology.nodes[face.n[2]].p;
        face.centroid = calculateFaceCentroid(p0, p1, p2).normalize();
    }
}

function reorderTriangleNodes(topology){

    for (var i = 0; i < topology.nodes.length; ++i)
    {
        var node = topology.nodes[i];
        var faceIndex = node.f[0];
        for (var j = 1; j < node.f.length - 1; ++j)
        {
            faceIndex = findNextFaceIndex(topology, i, faceIndex);
            var k = node.f.indexOf(faceIndex);
            node.f[k] = node.f[j];
            node.f[j] = faceIndex;
        }
    }
}

function distortMesh(mesh, degree){
    var totalSurfaceArea = 4 * Math.PI;
    var idealFaceArea = totalSurfaceArea / mesh.faces.length;
    var idealEdgeLength = Math.sqrt(idealFaceArea * 4 / Math.sqrt(3));
    var idealFaceHeight = idealEdgeLength * Math.sqrt(3) / 2;

    var rotationPredicate = function(oldNode0, oldNode1, newNode0, newNode1)
    {
        if (newNode0.f.length >= 7 ||
            newNode1.f.length >= 7 ||
            oldNode0.f.length <= 5 ||
            oldNode1.f.length <= 5) return false;
        var oldEdgeLength = oldNode0.p.distanceTo(oldNode1.p);
        var newEdgeLength = newNode0.p.distanceTo(newNode1.p);
        var ratio = oldEdgeLength / newEdgeLength;
        if (ratio >= 2 || ratio <= 0.5) return false;
        var v0 = oldNode1.p.clone().sub(oldNode0.p).divideScalar(oldEdgeLength);
        var v1 = newNode0.p.clone().sub(oldNode0.p).normalize();
        var v2 = newNode1.p.clone().sub(oldNode0.p).normalize();
        if (v0.dot(v1) < 0.2 || v0.dot(v2) < 0.2) return false;
        v0.negate();
        var v3 = newNode0.p.clone().sub(oldNode1.p).normalize();
        var v4 = newNode1.p.clone().sub(oldNode1.p).normalize();
        if (v0.dot(v3) < 0.2 || v0.dot(v4) < 0.2) return false;
        return true;
    };
    
    var i = 0;
    for(var i = 0;i<degree;i++){
        if (i >= degree) return;
        
        var consecutiveFailedAttempts = 0;
        var edgeIndex = Math.floor(Math.random()*mesh.edges.length);
        while (!conditionalRotateEdge(mesh, edgeIndex, rotationPredicate))
        {
            if (++consecutiveFailedAttempts >= mesh.edges.length) return false;
            edgeIndex = (edgeIndex + 1) % mesh.edges.length;
        }
    }

    return true;
}

function conditionalRotateEdge(mesh, edgeIndex, predicate){

    //var edgeIndex = Math.floor(Math.random()*planet.topology.edges.length);
    var edge = mesh.edges[edgeIndex];
    var face0 = mesh.faces[edge.f[0]];
    var face1 = mesh.faces[edge.f[1]];
    var farNodeFaceIndex0 = getFaceOppositeNodeIndex(face0, edge);
    var farNodeFaceIndex1 = getFaceOppositeNodeIndex(face1, edge);
    var newNodeIndex0 = face0.n[farNodeFaceIndex0];
    var oldNodeIndex0 = face0.n[(farNodeFaceIndex0 + 1) % 3];
    var newNodeIndex1 = face1.n[farNodeFaceIndex1];
    var oldNodeIndex1 = face1.n[(farNodeFaceIndex1 + 1) % 3];
    var oldNode0 = mesh.nodes[oldNodeIndex0];
    var oldNode1 = mesh.nodes[oldNodeIndex1];
    var newNode0 = mesh.nodes[newNodeIndex0];
    var newNode1 = mesh.nodes[newNodeIndex1];
    var newEdgeIndex0 = face1.e[(farNodeFaceIndex1 + 2) % 3];
    var newEdgeIndex1 = face0.e[(farNodeFaceIndex0 + 2) % 3];
    var newEdge0 = mesh.edges[newEdgeIndex0];
    var newEdge1 = mesh.edges[newEdgeIndex1];

    if (!predicate(oldNode0, oldNode1, newNode0, newNode1)) return false;
    
    oldNode0.e.splice(oldNode0.e.indexOf(edgeIndex), 1);
    oldNode1.e.splice(oldNode1.e.indexOf(edgeIndex), 1);
    newNode0.e.push(edgeIndex);
    newNode1.e.push(edgeIndex);

    edge.n[0] = newNodeIndex0;
    edge.n[1] = newNodeIndex1;
    
    newEdge0.f.splice(newEdge0.f.indexOf(edge.f[1]), 1);
    newEdge1.f.splice(newEdge1.f.indexOf(edge.f[0]), 1);
    newEdge0.f.push(edge.f[0]);
    newEdge1.f.push(edge.f[1]);
    
    oldNode0.f.splice(oldNode0.f.indexOf(edge.f[1]), 1);
    oldNode1.f.splice(oldNode1.f.indexOf(edge.f[0]), 1);
    newNode0.f.push(edge.f[1]);
    newNode1.f.push(edge.f[0]);
    
    face0.n[(farNodeFaceIndex0 + 2) % 3] = newNodeIndex1;
    face1.n[(farNodeFaceIndex1 + 2) % 3] = newNodeIndex0;

    face0.e[(farNodeFaceIndex0 + 1) % 3] = newEdgeIndex0;
    face1.e[(farNodeFaceIndex1 + 1) % 3] = newEdgeIndex1;
    face0.e[(farNodeFaceIndex0 + 2) % 3] = edgeIndex;
    face1.e[(farNodeFaceIndex1 + 2) % 3] = edgeIndex;

    return true;
}

function relaxMesh(mesh, multiplier){
    var totalSurfaceArea = 4 * Math.PI;
    var idealFaceArea = totalSurfaceArea / mesh.faces.length;
    var idealEdgeLength = Math.sqrt(idealFaceArea * 4 / Math.sqrt(3));
    var idealDistanceToCentroid = idealEdgeLength * Math.sqrt(3) / 3 * 0.9;
    
    var pointShifts = new Array(mesh.nodes.length);
    for (var i = 0; i < mesh.nodes.length; ++i)
        pointShifts[i] = new Vector3(0, 0, 0);

    for (var i = 0; i < mesh.faces.length; ++i){

        var face = mesh.faces[i];
        var n0 = mesh.nodes[face.n[0]];
        var n1 = mesh.nodes[face.n[1]];
        var n2 = mesh.nodes[face.n[2]];
        var p0 = n0.p;
        var p1 = n1.p;
        var p2 = n2.p;
        var e0 = p1.distanceTo(p0) / idealEdgeLength;
        var e1 = p2.distanceTo(p1) / idealEdgeLength;
        var e2 = p0.distanceTo(p2) / idealEdgeLength;
        var centroid = calculateFaceCentroid(p0, p1, p2).normalize();
        var v0 = centroid.clone().sub(p0);
        var v1 = centroid.clone().sub(p1);
        var v2 = centroid.clone().sub(p2);
        var length0 = v0.length();
        var length1 = v1.length();
        var length2 = v2.length();
        v0.multiplyScalar(multiplier * (length0 - idealDistanceToCentroid) / length0);
        v1.multiplyScalar(multiplier * (length1 - idealDistanceToCentroid) / length1);
        v2.multiplyScalar(multiplier * (length2 - idealDistanceToCentroid) / length2);
        pointShifts[face.n[0]].add(v0);
        pointShifts[face.n[1]].add(v1);
        pointShifts[face.n[2]].add(v2);
    }
    
    
    var origin = new Vector3(0, 0, 0);
    var plane = new THREE.Plane();

        for (var i = 0; i < mesh.nodes.length; ++i)
        {
            plane.setFromNormalAndCoplanarPoint(mesh.nodes[i].p, origin);
            pointShifts[i] = mesh.nodes[i].p.clone().add(plane.projectPoint(pointShifts[i])).normalize();
        }

    var rotationSupressions = new Array(mesh.nodes.length);
    for (var i = 0; i < mesh.nodes.length; ++i)
        rotationSupressions[i] = 0;
    
    var i = 0;
    for (var i = 0; i < mesh.edges.length; ++i){
        var edge = mesh.edges[i];
        var oldPoint0 = mesh.nodes[edge.n[0]].p;
        var oldPoint1 = mesh.nodes[edge.n[1]].p;
        var newPoint0 = pointShifts[edge.n[0]];
        var newPoint1 = pointShifts[edge.n[1]];
        var oldVector = oldPoint1.clone().sub(oldPoint0).normalize();
        var newVector = newPoint1.clone().sub(newPoint0).normalize();
        var suppression = (1 - oldVector.dot(newVector)) * 0.5;
        rotationSupressions[edge.n[0]] = Math.max(rotationSupressions[edge.n[0]], suppression);
        rotationSupressions[edge.n[1]] = Math.max(rotationSupressions[edge.n[1]], suppression);
    }

    
    var totalShift = 0;

        for (var i = 0; i < mesh.nodes.length; ++i)
        {
            var node = mesh.nodes[i];
            var point = node.p;
            var delta = point.clone();
            point.lerp(pointShifts[i], 1 - Math.sqrt(rotationSupressions[i])).normalize();
            delta.sub(point);
            totalShift += delta.length();
        }

    
    return totalShift;
}

function generatePlanetTiledTopology(topology){
    var corners = new Array(topology.faces.length);
    var borders = new Array(topology.edges.length);
    var tiles = new Array(topology.nodes.length);
    

    for (var i = 0; i < topology.faces.length; ++i){
        var face = topology.faces[i];
        corners[i] = new Corner(i, face.centroid.clone().multiplyScalar(1000), face.e.length, face.e.length, face.n.length);
    }

    for (var i = 0; i < topology.edges.length; ++i){
        var edge = topology.edges[i];
        borders[i] = new Border(i, 2, 4, 2); //edge.f.length, topology.faces[edge.f[0]].e.length + topology.faces[edge.f[1]].e.length - 2, edge.n.length
    }

    for (var i = 0; i < topology.nodes.length; ++i){
        var node = topology.nodes[i];
        tiles[i] = new Tile(i, node.p.clone().multiplyScalar(1000), node.f.length, node.e.length, node.e.length);
    }

    for (var i = 0; i < corners.length; ++i){
        var corner = corners[i];
        var face = topology.faces[i];
        for (var j = 0; j < face.e.length; ++j){
            corner.borders[j] = borders[face.e[j]];
        }
        for (var j = 0; j < face.n.length; ++j){
             corner.tiles[j] = tiles[face.n[j]];
        }
    }

    for (var i = 0; i < borders.length; ++i){
        var border = borders[i];
        var edge = topology.edges[i];
        var averageCorner = new Vector3(0, 0, 0);
        var n = 0;
        for (var j = 0; j < edge.f.length; ++j){
            var corner = corners[edge.f[j]]
            averageCorner.add(corner.position);
            border.corners[j] = corner;
            for (var k = 0; k < corner.borders.length; ++k){
                if (corner.borders[k] !== border) border.borders[n++] = corner.borders[k];
            }
        }
        border.midpoint = averageCorner.multiplyScalar(1 / border.corners.length);
        for (var j = 0; j < edge.n.length; ++j){
            border.tiles[j] = tiles[edge.n[j]];
        }
    }

    for (var i = 0; i < corners.length; ++i){
        var corner = corners[i];
        for (var j = 0; j < corner.borders.length; ++j){
            corner.corners[j] = corner.borders[j].oppositeCorner(corner);
        }
    }

    for (var i = 0; i < tiles.length; ++i){
        var tile = tiles[i];
        var node = topology.nodes[i];
        for (var j = 0; j < node.f.length; ++j){
            tile.corners[j] = corners[node.f[j]];
        }
        for (var j = 0; j < node.e.length; ++j){
            var border = borders[node.e[j]];
            if (border.tiles[0] === tile){
                for (var k = 0; k < tile.corners.length; ++k){
                    var corner0 = tile.corners[k];
                    var corner1 = tile.corners[(k + 1) % tile.corners.length];
                    if (border.corners[1] === corner0 && border.corners[0] === corner1){
                        border.corners[0] = corner0;
                        border.corners[1] = corner1;
                    }
                    else if (border.corners[0] !== corner0 || border.corners[1] !== corner1){
                        continue;
                    }
                    tile.borders[k] = border;
                    tile.tiles[k] = border.oppositeTile(tile);
                    break;
                }
            }
            else{
                for (var k = 0; k < tile.corners.length; ++k){
                    var corner0 = tile.corners[k];
                    var corner1 = tile.corners[(k + 1) % tile.corners.length];
                    if (border.corners[0] === corner0 && border.corners[1] === corner1){
                        border.corners[1] = corner0;
                        border.corners[0] = corner1;
                    }
                    else if (border.corners[1] !== corner0 || border.corners[0] !== corner1){
                        continue;
                    }
                    tile.borders[k] = border;
                    tile.tiles[k] = border.oppositeTile(tile);
                    break;
                }
            }
        }

        tile.averagePosition = new Vector3(0, 0, 0);
        for (var j = 0; j < tile.corners.length; ++j){
            tile.averagePosition.add(tile.corners[j].position);
        }
        tile.averagePosition.multiplyScalar(1 / tile.corners.length);
        
        var maxDistanceToCorner = 0;
        for (var j = 0; j < tile.corners.length; ++j){
            maxDistanceToCorner = Math.max(maxDistanceToCorner, tile.corners[j].position.distanceTo(tile.averagePosition));
        }
        
        var area = 0;
        for (var j = 0; j < tile.borders.length; ++j){
            area += calculateTriangleArea(tile.position, tile.borders[j].corners[0].position, tile.borders[j].corners[1].position);
        }

        tile.area = area;
        
        tile.normal = tile.position.clone().normalize();
        
        tile.boundingSphere = new THREE.Sphere(tile.averagePosition, maxDistanceToCorner);
    }

    for (var i = 0; i < corners.length; ++i){
        var corner = corners[i];
        corner.area = 0;
        for (var j = 0; j < corner.tiles.length; ++j){
            corner.area += corner.tiles[j].area / corner.tiles[j].corners.length;
        }
    }

    
   return  { corners: corners, borders: borders, tiles: tiles };
}



//TERRAIN
function generatePlanetTerrain(planet, plateCount,oceanicRate){

    var elevationBorderQueueSorter = function(left, right) { return left.distanceToPlateBoundary - right.distanceToPlateBoundary; };
    tiles = planet.tiledTopology.tiles;

    var plates = generateTectonicPlates(planet,planet.tiledTopology,plateCount,oceanicRate)
    identifyBoundaryBorders(planet.tiledTopology.borders)
    var corners = collectBoundaryCorners(planet.tiledTopology.corners)
    var indexes = calculatePlateBoundaryStress(corners)


    blurPlateBoundaryStress(corners,3,0.4)

    var elevationQueue = populateElevationBorderQueue(corners, indexes)
    processElevationBorderQueue(elevationQueue, elevationBorderQueueSorter);

    calculateTileAverageElevations(planet.tiledTopology.tiles);
}

function generateTectonicPlates(planet,topology, plateCount, oceanicRate){
    var plates = [];
    var platelessTiles = [];
    var platelessTilePlates = [];

    var failedCount = 0;
    while (plates.length < plateCount && failedCount < 10000){
        var corner = topology.corners[Math.floor(Math.random()*topology.corners.length)]
        var adjacentToExistingPlate = false;
        for (var i = 0; i < corner.tiles.length; ++i){
            if (corner.tiles[i].plate )
            {
                adjacentToExistingPlate = true;
                failedCount += 1;
                break;
            }
        }
        if (adjacentToExistingPlate) continue;
        
        failedCount = 0;
        
        var oceanic = (Math.random() < oceanicRate);
        var plate = new Plate(
                new THREE.Color("#"+((1<<24)*Math.random()|0).toString(16)),
                new randomUnitVector(Math.random()),
                randomFromInterval(-Math.PI / 30, Math.PI / 30),
                randomFromInterval(-Math.PI / 30, Math.PI / 30),
                oceanic ? -randomFromInterval(0.8, 0.3) : randomFromInterval(0.1, 0.5),
                oceanic,
                corner);
            
        plates.push(plate);

        for (var i = 0; i < corner.tiles.length; ++i){
            corner.tiles[i].plate = plate;
            plate.tiles.push(corner.tiles[i]);
        }

        for (var i = 0; i < corner.tiles.length; ++i){
            var corner = corner.tiles[i];
            for (var j = 0; j < corner.tiles.length; ++j)
            {
                var adjacentTile = corner.tiles[j];
                if (!adjacentTile.plate)
                {
                    platelessTiles.push(adjacentTile);
                    platelessTilePlates.push(plate);
                }
            }
        }
    }
    while (platelessTiles.length > 0){
        var tileIndex = Math.floor(Math.pow(Math.random(), 2) * platelessTiles.length);
        var tile = platelessTiles[tileIndex];
        var plate = platelessTilePlates[tileIndex];
        platelessTiles.splice(tileIndex, 1);
        platelessTilePlates.splice(tileIndex, 1);
        if (!tile.plate)
        {
            tile.plate = plate;
            plate.tiles.push(tile);
            for (var j = 0; j < tile.tiles.length; ++j)
            {
                if (!tile.tiles[j].plate)
                {
                    platelessTiles.push(tile.tiles[j]);
                    platelessTilePlates.push(plate);
                }
            }
        }
    }
    
    calculateCornerDistancesToPlateRoot(plates)
    planet.plates = plates
    return plates
}

function calculateCornerDistancesToPlateRoot(plates){
    var distanceCornerQueue = [];
    for (var i = 0; i < plates.length; ++i)
    {
        var corner = plates[i].root;
        corner.distanceToPlateRoot = 0;
        for (var j = 0; j < corner.corners.length; ++j)
        {
            distanceCornerQueue.push({ corner: corner.corners[j], distanceToPlateRoot: corner.borders[j].length() });
        }
    }
    
    var distanceCornerQueueSorter = function(left, right) { return left.distanceToPlateRoot - right.distanceToPlateRoot; };

    while(true){
        if (distanceCornerQueue.length === 0) return;

        var iEnd = iEnd = distanceCornerQueue.length;
        for (var i = 0; i < iEnd; ++i)
        {
            var front = distanceCornerQueue[i];
            var corner = front.corner;
            var distanceToPlateRoot = front.distanceToPlateRoot;
            if (!corner.distanceToPlateRoot || corner.distanceToPlateRoot > distanceToPlateRoot)
            {
                corner.distanceToPlateRoot = distanceToPlateRoot;
                for (var j = 0; j < corner.corners.length; ++j)
                {
                    distanceCornerQueue.push({ corner: corner.corners[j], distanceToPlateRoot: distanceToPlateRoot + corner.borders[j].length() });
                }
            }
        }
        distanceCornerQueue.splice(0, iEnd);
        distanceCornerQueue.sort(distanceCornerQueueSorter);  
    };
}

function identifyBoundaryBorders(borders){
    for (var i = 0; i < borders.length; ++i)
    {
        var border = borders[i];
        if (border.tiles[0].plate !== border.tiles[1].plate)
        {
            border.betweenPlates = true;
            border.corners[0].betweenPlates = true;
            border.corners[1].betweenPlates = true;
            border.tiles[0].plate.boundaryBorders.push(border);
            border.tiles[1].plate.boundaryBorders.push(border);
        }
    }
}

function drawBoundaryBorders(plates){
    for (var i =0;i<plates.length;i++){
        var colorDeviance = new THREE.Color(0.2,0.2,0.2);
        var material = new THREE.LineBasicMaterial({});
        var geometry = new THREE.Geometry();
        var line = new THREE.Line(geometry, material);

        plates[i].boundaryBorders.forEach(function(x){
            var normal = x.midpoint.clone().normalize();
            var offset = normal.clone().multiplyScalar(1);

            var a = x.corners[0].position;
            var b = x.corners[1].position;

            var p = x.corners[0].pressure;
            geometry.vertices.push(a.clone().add(offset));
            geometry.vertices.push(b.clone().add(offset)); 

            var pressure = Math.max(-1, Math.min((x.corners[0].pressure + x.corners[1].pressure) / 2, 1));
            material.color =  (pressure <= 0) ? new THREE.Color(1 + pressure, 1, 0) : new THREE.Color(1, 1 - pressure, 0);
            var line = new THREE.Line(geometry, material);
            scene.add(line);
            geometry = new THREE.Geometry();
        })  
    }
}

function collectBoundaryCorners(corners){
    var boundaryCorners = [];
    for (var j = 0; j < corners.length; ++j)
    {
        var corner = corners[j];
        if (corner.betweenPlates)
        {
            boundaryCorners.push(corner);
            corner.tiles[0].plate.boundaryCorners.push(corner);
            if (corner.tiles[1].plate !== corner.tiles[0].plate) corner.tiles[1].plate.boundaryCorners.push(corner);
            if (corner.tiles[2].plate !== corner.tiles[0].plate && corner.tiles[2].plate !== corner.tiles[1].plate) corner.tiles[2].plate.boundaryCorners.push(corner);
        }
    }
    
    return boundaryCorners;
}

function calculatePlateBoundaryStress(boundaryCorners){
    var boundaryCornerInnerBorderIndexes = new Array(boundaryCorners.length);
    for (var i = 0; i < boundaryCorners.length; ++i)
    {
        var corner = boundaryCorners[i];
        corner.distanceToPlateBoundary = 0;
    
        var innerBorder;
        var innerBorderIndex;
        for (var j = 0; j < corner.borders.length; ++j){
            var border = corner.borders[j];
            if (!border.betweenPlates){
                innerBorder = border;
                innerBorderIndex = j;
                break;
            }
        }
        
        if (innerBorder){
            boundaryCornerInnerBorderIndexes[i] = innerBorderIndex;
            var outerBorder0 = corner.borders[(innerBorderIndex + 1) % corner.borders.length];
            var outerBorder1 = corner.borders[(innerBorderIndex + 2) % corner.borders.length]
            var farCorner0 = outerBorder0.oppositeCorner(corner);
            var farCorner1 = outerBorder1.oppositeCorner(corner);
            var plate0 = innerBorder.tiles[0].plate;
            var plate1 = outerBorder0.tiles[0].plate !== plate0 ? outerBorder0.tiles[0].plate : outerBorder0.tiles[1].plate;
            var boundaryVector = farCorner0.vectorTo(farCorner1);
            var boundaryNormal = boundaryVector.clone().cross(corner.position);
            var stress = calculateStress(plate0.calculateMovement(corner.position), plate1.calculateMovement(corner.position), boundaryVector, boundaryNormal);
            corner.pressure = stress.pressure;
            corner.shear = stress.shear;
        }
        else{
            boundaryCornerInnerBorderIndexes[i] = null;
            var plate0 = corner.tiles[0].plate;
            var plate1 = corner.tiles[1].plate;
            var plate2 = corner.tiles[2].plate;
            var boundaryVector0 = corner.corners[0].vectorTo(corner);
            var boundaryVector1 = corner.corners[1].vectorTo(corner);
            var boundaryVector2 = corner.corners[2].vectorTo(corner);
            var boundaryNormal0 = boundaryVector0.clone().cross(corner.position);
            var boundaryNormal1 = boundaryVector1.clone().cross(corner.position);
            var boundaryNormal2 = boundaryVector2.clone().cross(corner.position);
            var stress0 = calculateStress(plate0.calculateMovement(corner.position), plate1.calculateMovement(corner.position), boundaryVector0, boundaryNormal0);
            var stress1 = calculateStress(plate1.calculateMovement(corner.position), plate2.calculateMovement(corner.position), boundaryVector1, boundaryNormal1);
            var stress2 = calculateStress(plate2.calculateMovement(corner.position), plate0.calculateMovement(corner.position), boundaryVector2, boundaryNormal2);
            
            corner.pressure = (stress0.pressure + stress1.pressure + stress2.pressure) / 3;
            corner.shear = (stress0.shear + stress1.shear + stress2.shear) / 3;
        }
    }
    
    return boundaryCornerInnerBorderIndexes;
}

function calculateStress(movement0, movement1, boundaryVector, boundaryNormal){
    var relativeMovement = movement0.clone().sub(movement1);
    var pressureVector = relativeMovement.clone().projectOnVector(boundaryNormal);
    var pressure = pressureVector.length();
    if (pressureVector.dot(boundaryNormal) > 0) pressure = -pressure;
    var shear = relativeMovement.clone().projectOnVector(boundaryVector).length();
    return { pressure: 2 / (1 + Math.exp(-pressure / 30)) - 1, shear: 2 / (1 + Math.exp(-shear / 30)) - 1 };
}

function blurPlateBoundaryStress(boundaryCorners, stressBlurIterations, stressBlurCenterWeighting, action){
    var newCornerPressure = new Array(boundaryCorners.length);
    var newCornerShear = new Array(boundaryCorners.length);
    for (var i = 0; i < stressBlurIterations; ++i)
    {
        for (var j = 0; j < boundaryCorners.length; ++j)
        {
            var corner = boundaryCorners[j];
            var averagePressure = 0;
            var averageShear = 0;
            var neighborCount = 0;
            for (var k = 0; k < corner.corners.length; ++k)
            {
                var neighbor = corner.corners[k];
                if (neighbor.betweenPlates)
                {
                    averagePressure += neighbor.pressure;
                    averageShear += neighbor.shear;
                    ++neighborCount;
                }
            }
            newCornerPressure[j] = corner.pressure * stressBlurCenterWeighting + (averagePressure / neighborCount) * (1 - stressBlurCenterWeighting);
            newCornerShear[j] = corner.shear * stressBlurCenterWeighting + (averageShear / neighborCount) * (1 - stressBlurCenterWeighting);
        }
        
        for (var j = 0; j < boundaryCorners.length; ++j)
        {
            var corner = boundaryCorners[j];
            if (corner.betweenPlates)
            {
                corner.pressure = newCornerPressure[j];
                corner.shear = newCornerShear[j];
            }
        }
    }
}

function populateElevationBorderQueue(boundaryCorners, boundaryCornerInnerBorderIndexes){
    var elevationBorderQueue = [];
    for (var i = 0; i < boundaryCorners.length; ++i)
    {
        var corner = boundaryCorners[i];
        
        var innerBorderIndex = boundaryCornerInnerBorderIndexes[i];
        if (innerBorderIndex !== null)
        {
            var innerBorder = corner.borders[innerBorderIndex];
            var outerBorder0 = corner.borders[(innerBorderIndex + 1) % corner.borders.length];
            var plate0 = innerBorder.tiles[0].plate;
            var plate1 = outerBorder0.tiles[0].plate !== plate0 ? outerBorder0.tiles[0].plate : outerBorder0.tiles[1].plate;
            
            var calculateElevation;
            
            if (corner.pressure > 0.3)
            {
                corner.elevation = Math.max(plate0.elevation, plate1.elevation) + corner.pressure;
                if (plate0.oceanic === plate1.oceanic)
                    calculateElevation = calculateCollidingElevation;
                else if (plate0.oceanic)
                    calculateElevation = calculateSubductingElevation;
                else
                    calculateElevation = calculateSuperductingElevation;
            }
            else if (corner.pressure < -0.3)
            {
                corner.elevation = Math.max(plate0.elevation, plate1.elevation) - corner.pressure / 4;
                calculateElevation = calculateDivergingElevation;
            }
            else if (corner.shear > 0.3)
            {
                corner.elevation = Math.max(plate0.elevation, plate1.elevation) + corner.shear / 8;
                calculateElevation = calculateShearingElevation;
            }
            else
            {
                corner.elevation = (plate0.elevation + plate1.elevation) / 2;
                calculateElevation = calculateDormantElevation;
            }
            
            var nextCorner = innerBorder.oppositeCorner(corner);
            if (!nextCorner.betweenPlates)
            {
                elevationBorderQueue.push({
                    origin: {
                        corner: corner,
                        pressure: corner.pressure,
                        shear: corner.shear,
                        plate: plate0,
                        calculateElevation: calculateElevation },
                    border: innerBorder,
                    corner: corner,
                    nextCorner: nextCorner,
                    distanceToPlateBoundary: innerBorder.length(),
                });
            }
        }
        else
        {
            var plate0 = corner.tiles[0].plate;
            var plate1 = corner.tiles[1].plate;
            var plate2 = corner.tiles[2].plate;
            
            elevation = 0;

            if (corner.pressure > 0.3)
            {
                corner.elevation = Math.max(plate0.elevation, plate1.elevation, plate2.elevation) + corner.pressure;
            }
            else if (corner.pressure < -0.3)
            {
                corner.elevation = Math.max(plate0.elevation, plate1.elevation, plate2.elevation) + corner.pressure / 4;
            }
            else if (corner.shear > 0.3)
            {
                corner.elevation = Math.max(plate0.elevation, plate1.elevation, plate2.elevation) + corner.shear / 8;
            }
            else
            {
                corner.elevation = (plate0.elevation + plate1.elevation + plate2.elevation) / 3;
            }
        }
    }
    
    return elevationBorderQueue;
}

function calculateCollidingElevation(distanceToPlateBoundary, distanceToPlateRoot, boundaryElevation, plateElevation, pressure, shear){
    var t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
    if (t < 0.5)
    {
        t = t / 0.5;
        return plateElevation + Math.pow(t - 1, 2) * (boundaryElevation - plateElevation);
    }
    else
    {
        return plateElevation;
    }
}

function calculateSuperductingElevation(distanceToPlateBoundary, distanceToPlateRoot, boundaryElevation, plateElevation, pressure, shear){
    var t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
    if (t < 0.2)
    {
        t = t / 0.2;
        return boundaryElevation + t * (plateElevation - boundaryElevation + pressure / 2);
    }
    else if (t < 0.5)
    {
        t = (t - 0.2) / 0.3;
        return plateElevation + Math.pow(t - 1, 2) * pressure / 2;
    }
    else
    {
        return plateElevation;
    }
}

function calculateSubductingElevation(distanceToPlateBoundary, distanceToPlateRoot, boundaryElevation, plateElevation, pressure, shear){
    var t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
    return plateElevation + Math.pow(t - 1, 2) * (boundaryElevation - plateElevation);
}

function calculateDivergingElevation(distanceToPlateBoundary, distanceToPlateRoot, boundaryElevation, plateElevation, pressure, shear){
    var t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
    if (t < 0.3)
    {
        t = t / 0.3;
        return plateElevation + Math.pow(t - 1, 2) * (boundaryElevation - plateElevation);
    }
    else
    {
        return plateElevation;
    }
}

function calculateShearingElevation(distanceToPlateBoundary, distanceToPlateRoot, boundaryElevation, plateElevation, pressure, shear){
    var t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
    if (t < 0.2)
    {
        t = t / 0.2;
        return plateElevation + Math.pow(t - 1, 2) * (boundaryElevation - plateElevation);
    }
    else
    {
        return plateElevation;
    }
}

function calculateDormantElevation(distanceToPlateBoundary, distanceToPlateRoot, boundaryElevation, plateElevation, pressure, shear){
    var t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
    var elevationDifference = boundaryElevation - plateElevation;
    var a = 2 * elevationDifference;
    var b = -3 * elevationDifference;
    return t * t * elevationDifference * (2 * t - 3) + boundaryElevation;
}

function processElevationBorderQueue(elevationBorderQueue, elevationBorderQueueSorter){
    while(true){
        if (elevationBorderQueue.length === 0) return;

        var iEnd = iEnd = elevationBorderQueue.length;
        for (var i = 0; i < iEnd; ++i)
        {
            var front = elevationBorderQueue[i];
            var corner = front.nextCorner;
            if (!corner.elevation)
            {
                corner.distanceToPlateBoundary = front.distanceToPlateBoundary;
                
                corner.elevation = front.origin.calculateElevation(
                    corner.distanceToPlateBoundary,
                    corner.distanceToPlateRoot,
                    front.origin.corner.elevation,
                    front.origin.plate.elevation,
                    front.origin.pressure,
                    front.origin.shear);
 
                for (var j = 0; j < corner.borders.length; ++j)
                {
                    var border = corner.borders[j];
                    if (!border.betweenPlates)
                    {
                        var nextCorner = corner.corners[j];
                        var distanceToPlateBoundary = corner.distanceToPlateBoundary + border.length();
                        if (!nextCorner.distanceToPlateBoundary || nextCorner.distanceToPlateBoundary > distanceToPlateBoundary)
                        {
                            elevationBorderQueue.push({
                                origin: front.origin,
                                border: border,
                                corner: corner,
                                nextCorner: nextCorner,
                                distanceToPlateBoundary: distanceToPlateBoundary,
                            });
                        }
                    }
                }
            }
        }
        elevationBorderQueue.splice(0, iEnd);
        elevationBorderQueue.sort(elevationBorderQueueSorter);
    }
}

function calculateTileAverageElevations(tiles){

    for (var i = 0; i < tiles.length; ++i){
        var tile = tiles[i];
        var elevation = 0;
        for (var j = 0; j < tile.corners.length; ++j){
            elevation += tile.corners[j].elevation;
        }
        tile.elevation = elevation / tile.corners.length;
    }
}



//WEATHER
function generatePlanetWeather(topology, heatLevel, moistureLevel){
    var planetRadius = 1000;
    var whorls;
    var activeCorners;
    var totalHeat;
    var remainingHeat;
    var totalMoisture;
    var remainingMoisture;

    whorls = generateAirCurrentWhorls(planetRadius)
    calculateAirCurrents(topology.corners, whorls, planetRadius)


    initialHeat = initializeAirHeat(topology.corners,heatLevel )

    activeCorners = initialHeat.corners;
    totalHeat = initialHeat.airHeat;
    remainingHeat = initialHeat.airHeat;

    var consumedHeat = processAirHeat(activeCorners);
    while(remainingHeat > 0 && consumedHeat >= 0.0001){
        var consumedHeat = processAirHeat(activeCorners);
        remainingHeat -= consumedHeat;
    }
    calculateTemperature(topology.corners, topology.tiles, planetRadius);


    initialMoisture = initializeAirMoisture(topology.corners, moistureLevel)
    

    activeCorners = initialMoisture.corners;
    totalMoisture = initialMoisture.airMoisture;
    remainingMoisture = initialMoisture.airMoisture;
    var consumedMoisture = processAirMoisture(activeCorners);

    while(remainingMoisture > 0 && consumedMoisture >= 0.0001){
        var consumedMoisture = processAirMoisture(activeCorners);
        remainingMoisture -= consumedMoisture;
   }
    calculateMoisture(topology.corners, topology.tiles)
}

function generateAirCurrentWhorls(planetRadius){
    var random = Math.random()
    var whorls = [];
    var direction = randomInteger(0, 1) ? 1 : -1;
    var layerCount = randomInteger(4, 7);
    var circumference = Math.PI * 2 * planetRadius;
    var fullRevolution = Math.PI * 2;
    var baseWhorlRadius = circumference / (2 * (layerCount - 1));

    whorls.push({
        center: new Vector3(0, planetRadius, 0)
            .applyAxisAngle(new Vector3(1, 0, 0), randomRealInclusive(0, fullRevolution / (2 * (layerCount + 4))))
            .applyAxisAngle(new Vector3(0, 1, 0), randomReal(0, fullRevolution)),
        strength: randomRealInclusive(fullRevolution / 36, fullRevolution / 24) * direction,
        radius: randomRealInclusive(baseWhorlRadius * 0.8, baseWhorlRadius * 1.2) });

    for (var i = 1; i < layerCount - 1; ++i)
    {
        direction = -direction;
        var baseTilt = i / (layerCount - 1) * fullRevolution / 2;
        var layerWhorlCount = Math.ceil((Math.sin(baseTilt) * planetRadius * fullRevolution) / baseWhorlRadius);
        for (var j = 0; j < layerWhorlCount; ++j)
        {
            whorls.push({
                center: new Vector3(0, planetRadius, 0)
                    .applyAxisAngle(new Vector3(1, 0, 0), randomRealInclusive(0, fullRevolution / (2 * (layerCount + 4))))
                    .applyAxisAngle(new Vector3(0, 1, 0), randomReal(0, fullRevolution))
                    .applyAxisAngle(new Vector3(1, 0, 0), baseTilt)
                    .applyAxisAngle(new Vector3(0, 1, 0), fullRevolution * (j + (i % 2) / 2) / layerWhorlCount),
                strength: randomRealInclusive(fullRevolution / 48, fullRevolution / 32) * direction,
                radius: randomRealInclusive(baseWhorlRadius * 0.8, baseWhorlRadius * 1.2) });
        }
    }

    direction = -direction;
    whorls.push({
        center: new Vector3(0, planetRadius, 0)
            .applyAxisAngle(new Vector3(1, 0, 0), randomRealInclusive(0, fullRevolution / (2 * (layerCount + 4))))
            .applyAxisAngle(new Vector3(0, 1, 0), randomReal(0, fullRevolution))
            .applyAxisAngle(new Vector3(1, 0, 0), fullRevolution / 2),
        strength: randomRealInclusive(fullRevolution / 36, fullRevolution / 24) * direction,
        radius: randomRealInclusive(baseWhorlRadius * 0.8, baseWhorlRadius * 1.2) });
        
    return whorls;
}

function calculateAirCurrents(corners, whorls, planetRadius){
    var i = 0;
    for(var i = 0;i<corners.length;i++)
    {
        if (i >= corners.length) return;

        var corner = corners[i];
        var airCurrent = new Vector3(0, 0, 0);
        var weight = 0;
        for (var j = 0; j < whorls.length; ++j)
        {
            var whorl = whorls[j];
            var angle = whorl.center.angleTo(corner.position);
            var distance = angle * planetRadius;
            if (distance < whorl.radius)
            {
                var normalizedDistance = distance / whorl.radius;
                var whorlWeight = 1 - normalizedDistance;
                var whorlStrength = planetRadius * whorl.strength * whorlWeight * normalizedDistance;
                var whorlCurrent = whorl.center.clone().cross(corner.position).setLength(whorlStrength);
                airCurrent.add(whorlCurrent);
                weight += whorlWeight;
            }
        }
        airCurrent.divideScalar(weight);
        corner.airCurrent = airCurrent;
        corner.airCurrentSpeed = airCurrent.length(); //kilometers per hour
        
        corner.airCurrentOutflows = new Array(corner.borders.length);
        var airCurrentDirection = airCurrent.clone().normalize();
        var outflowSum = 0;
        for (var j = 0; j < corner.corners.length; ++j)
        {
            var vector = corner.vectorTo(corner.corners[j]).normalize();
            var dot = vector.dot(airCurrentDirection);
            if (dot > 0)
            {
                corner.airCurrentOutflows[j] = dot;
                outflowSum += dot;
            }
            else
            {
                corner.airCurrentOutflows[j] = 0;
            }
        }

        if (outflowSum > 0)
        {
            for (var j = 0; j < corner.borders.length; ++j)
            {
                corner.airCurrentOutflows[j] /= outflowSum;
            }
        }
    };
}

function initializeAirHeat(corners, heatLevel){
    var activeCorners = [];
    var airHeat = 0;
    for (var i = 0; i < corners.length; ++i)
    {
        var corner = corners[i];
        corner.airHeat = corner.area * heatLevel;
        corner.newAirHeat = 0;
        corner.heat = 0;

        corner.heatAbsorption = 0.1 * corner.area / Math.max(0.1, Math.min(corner.airCurrentSpeed, 1));
        if (corner.elevation <= 0)
        {
            corner.maxHeat = corner.area;
        }
        else
        {
            corner.maxHeat = corner.area;
            corner.heatAbsorption *= 2;
        }
        
        activeCorners.push(corner);
        airHeat += corner.airHeat;
    }
    
    return{ corners: activeCorners, airHeat: airHeat };
}

function processAirHeat(activeCorners){
    var consumedHeat = 0;
    var activeCornerCount = activeCorners.length;
    for (var i = 0; i < activeCornerCount; ++i)
    {
        var corner = activeCorners[i];
        if (corner.airHeat === 0) continue;
        
        var heatChange = Math.max(0, Math.min(corner.airHeat, corner.heatAbsorption * (1 - corner.heat / corner.maxHeat)));
        corner.heat += heatChange;
        consumedHeat += heatChange;
        var heatLoss = corner.area * (corner.heat / corner.maxHeat) * 0.02;
        heatChange = Math.min(corner.airHeat, heatChange + heatLoss);
        
        var remainingCornerAirHeat = corner.airHeat - heatChange;
        corner.airHeat = 0;
        
        for (var j = 0; j < corner.corners.length; ++j)
        {
            var outflow = corner.airCurrentOutflows[j];
            if (outflow > 0)
            {
                corner.corners[j].newAirHeat += remainingCornerAirHeat * outflow;
                activeCorners.push(corner.corners[j]);
            }
        }
    }
    
    activeCorners.splice(0, activeCornerCount);
    
    for (var i = 0; i < activeCorners.length; ++i)
    {
        var corner = activeCorners[i];
        corner.airHeat = corner.newAirHeat;
    }
    for (var i = 0; i < activeCorners.length; ++i)
    {
        activeCorners[i].newAirHeat = 0;
    }
    
    return consumedHeat;
}

function calculateTemperature(corners, tiles, planetRadius){
    for (var i = 0; i < corners.length; ++i)
    {
        var corner = corners[i];
        var latitudeEffect = Math.sqrt(1 - Math.abs(corner.position.y) / planetRadius );
        var elevationEffect = 1 - Math.pow(Math.max(0, Math.min(corner.elevation * 0.8, 1)), 2);
        var normalizedHeat = corner.heat / corner.area;
        corner.temperature = (latitudeEffect * elevationEffect * 0.7 + normalizedHeat * 0.3) * 5/3 - 2/3;
        
        delete corner.airHeat;
        delete corner.newAirHeat;
        delete corner.heat;
        delete corner.maxHeat;
        delete corner.heatAbsorption;
    }

    for (var i = 0; i < tiles.length; ++i)
    {
        var tile = tiles[i];
        tile.temperature = 0;
        for (var j = 0; j < tile.corners.length; ++j)
        {
            tile.temperature += tile.corners[j].temperature;
        }
        tile.temperature /= tile.corners.length;
    }
}

function initializeAirMoisture(corners, moistureLevel){
    activeCorners = [];
    airMoisture = 0;
    for (var i = 0; i < corners.length; ++i)
    {
        var corner = corners[i];
        corner.airMoisture = (corner.elevation > 0) ? 0 : corner.area * moistureLevel * Math.max(0, Math.min(0.5 + corner.temperature * 0.5, 1));
        corner.newAirMoisture = 0;
        corner.precipitation = 0;

        corner.precipitationRate = 0.0075 * corner.area / Math.max(0.1, Math.min(corner.airCurrentSpeed, 1));
        corner.precipitationRate *= 1 + (1 - Math.max(0, Math.max(corner.temperature, 1))) * 0.1;
        if (corner.elevation > 0)
        {
            corner.precipitationRate *= 1 + corner.elevation * 0.5;
            corner.maxPrecipitation = corner.area * (0.25 + Math.max(0, Math.min(corner.elevation, 1)) * 0.25);
        }
        else
        {
            corner.maxPrecipitation = corner.area * 0.25;
        }
        
        activeCorners.push(corner);
        airMoisture += corner.airMoisture;
    }
    return{ corners: activeCorners, airMoisture: airMoisture };
}

function processAirMoisture(activeCorners){
    var consumedMoisture = 0;
    var activeCornerCount = activeCorners.length;

    for (var i = 0; i < activeCornerCount; ++i)
    {
        var corner = activeCorners[i];

        if (corner.airMoisture === 0) continue;
        
        var moistureChange = Math.max(0, Math.min(corner.airMoisture, corner.precipitationRate * (1 - corner.precipitation / corner.maxPrecipitation)));
        corner.precipitation += moistureChange;
        consumedMoisture += moistureChange;
        var moistureLoss = corner.area * (corner.precipitation / corner.maxPrecipitation) * 0.02;
        moistureChange = Math.min(corner.airMoisture, moistureChange + moistureLoss);
        
        var remainingCornerAirMoisture = corner.airMoisture - moistureChange;
        corner.airMoisture = 0;
        
        for (var j = 0; j < corner.corners.length; ++j)
        {
            var outflow = corner.airCurrentOutflows[j];
            if (outflow > 0)
            {
                corner.corners[j].newAirMoisture += remainingCornerAirMoisture * outflow;
                activeCorners.push(corner.corners[j]);

            }
        }
    }

    activeCorners.splice(0, activeCornerCount);

    for (var i = 0; i < activeCorners.length; ++i)
    {
        var corner = activeCorners[i];
        corner.airMoisture = corner.newAirMoisture;
    }
    for (var i = 0; i < activeCorners.length; ++i)
    {
        activeCorners[i].newAirMoisture = 0;
    }
    
    return consumedMoisture;
}

function calculateMoisture(corners, tiles){
    for (var i = 0; i < corners.length; ++i)
    {
        var corner = corners[i];
        corner.moisture = corner.precipitation / corner.area / 0.5;
        delete corner.airMoisture;
        delete corner.newAirMoisture;
        delete corner.precipitation;
        delete corner.maxPrecipitation;
        delete corner.precipitationRate;
    }

    for (var i = 0; i < tiles.length; ++i)
    {
        var tile = tiles[i];
        tile.moisture = 0;
        for (var j = 0; j < tile.corners.length; ++j)
        {
            tile.moisture += tile.corners[j].moisture;
        }
        tile.moisture /= tile.corners.length;
    }
}

function generatePlanetBiomes(tiles, planetRadius){
    for (var i = 0; i < tiles.length; ++i)
    {
        var tile = tiles[i];
        var elevation = Math.max(0, tile.elevation);
        var latitude = Math.abs(tile.position.y / planetRadius);
        var temperature = tile.temperature;
        var moisture = tile.moisture;
        
        if (elevation <= 0)
        {
            if (temperature > 0)
            {
                tile.biome = "ocean";
            }
            else
            {
                tile.biome = "oceanGlacier";
            }
        }
        else if (elevation < 0.6)
        {
            if (temperature > 0.75)
            {
                if (moisture < 0.25)
                {
                    tile.biome = "desert";
                }
                else
                {
                    tile.biome = "rainForest";
                }
            }
            else if (temperature > 0.5)
            {
                if (moisture < 0.25)
                {
                    tile.biome = "rocky";
                }
                else if (moisture < 0.50)
                {
                    tile.biome = "plains";
                }
                else
                {
                    tile.biome = "swamp";
                }
            }
            else if (temperature > 0)
            {
                if (moisture < 0.25)
                {
                    tile.biome = "plains";
                }
                else if (moisture < 0.50)
                {
                    tile.biome = "grassland";
                }
                else
                {
                    tile.biome = "deciduousForest";
                }
            }
            else
            {
                if (moisture < 0.25)
                {
                    tile.biome = "tundra";
                }
                else
                {
                    tile.biome = "landGlacier";
                }
            }
        }
        else if (elevation < 0.8)
        {
            if (temperature > 0)
            {
                if (moisture < 0.25)
                {
                    tile.biome = "tundra";
                }
                else
                {
                    tile.biome = "coniferForest";
                }
            }
            else
            {
                tile.biome = "tundra";
            }
        }
        else
        {
            if (temperature > 0 || moisture < 0.25)
            {
                tile.biome = "mountain";
            }
            else
            {
                tile.biome = "snowyMountain";
            }
        }
    }
}



//RENDER
function generatePlanetRenderData(topology){
    var renderData = {};
    console.log("Generating Surface Visuals")
    renderData.surface = buildSurfaceRenderObject(topology)
    console.log("Generatingd Plate Boundaries Visuals")
    renderData.plateBoundaries = buildPlateBoundariesRenderObject(topology.borders)
    console.log("Generating Plate movement Visuals")
    renderData.plateMovements = buildPlateMovementsRenderObject(topology.tiles)
    console.log("Generating Air Current Visuals")
    renderData.airCurrents = buildAirCurrentsRenderObject(topology.corners)
   
    return renderData
}

function buildSurfaceRenderObject(topology){
    var planetGeometry = new THREE.Geometry();
    var terrainColors = [];
    var plateColors = [];
    var elevationColors = [];
    var temperatureColors = [];
    var moistureColors = [];
    
    var i = 0;
    for(var i = 0;i<topology.tiles.length;i++){

        var tile = tiles[i];
        
        var colorDeviance = new THREE.Color(randomUnit(), randomUnit(), randomUnit());
        var terrainColor;
        if (tile.elevation <= 0)
        {
            var normalizedElevation = Math.min(-tile.elevation, 1);
            if (tile.biome === "ocean") terrainColor = new THREE.Color(0x0066FF).lerp(new THREE.Color(0x0044BB), Math.min(-tile.elevation, 1)).lerp(colorDeviance, 0.10);
            else if (tile.biome === "oceanGlacier") terrainColor = new THREE.Color(0xDDEEFF).lerp(colorDeviance, 0.10);
            else terrainColor = new THREE.Color(0xFF00FF);
        }
        else if (tile.elevation < 0.6)
        {
            var normalizedElevation = tile.elevation / 0.6;
            if (tile.biome === "desert") terrainColor = new THREE.Color(0xDDDD77).lerp(new THREE.Color(0xBBBB55), normalizedElevation).lerp(colorDeviance, 0.10);
            else if (tile.biome === "rainForest") terrainColor = new THREE.Color(0x44DD00).lerp(new THREE.Color(0x229900), normalizedElevation).lerp(colorDeviance, 0.20);
            else if (tile.biome === "rocky") terrainColor = new THREE.Color(0xAA9977).lerp(new THREE.Color(0x887755), normalizedElevation).lerp(colorDeviance, 0.15);
            else if (tile.biome === "plains") terrainColor = new THREE.Color(0x99BB44).lerp(new THREE.Color(0x667722), normalizedElevation).lerp(colorDeviance, 0.10);
            else if (tile.biome === "grassland") terrainColor = new THREE.Color(0x77CC44).lerp(new THREE.Color(0x448822), normalizedElevation).lerp(colorDeviance, 0.15);
            else if (tile.biome === "swamp") terrainColor = new THREE.Color(0x77AA44).lerp(new THREE.Color(0x446622), normalizedElevation).lerp(colorDeviance, 0.25);
            else if (tile.biome === "deciduousForest") terrainColor = new THREE.Color(0x33AA22).lerp(new THREE.Color(0x116600), normalizedElevation).lerp(colorDeviance, 0.10);
            else if (tile.biome === "tundra") terrainColor = new THREE.Color(0x9999AA).lerp(new THREE.Color(0x777788), normalizedElevation).lerp(colorDeviance, 0.15);
            else if (tile.biome === "landGlacier") terrainColor = new THREE.Color(0xDDEEFF).lerp(colorDeviance, 0.10);
            else terrainColor = new THREE.Color(0xFF00FF);
        }
        else if (tile.elevation < 0.8)
        {
            var normalizedElevation = (tile.elevation - 0.6) / 0.2;
            if (tile.biome === "tundra") terrainColor = new THREE.Color(0x777788).lerp(new THREE.Color(0x666677), normalizedElevation).lerp(colorDeviance, 0.10);
            else if (tile.biome === "coniferForest") terrainColor = new THREE.Color(0x338822).lerp(new THREE.Color(0x116600), normalizedElevation).lerp(colorDeviance, 0.10);
            else if (tile.biome === "snow") terrainColor = new THREE.Color(0xEEEEEE).lerp(new THREE.Color(0xDDDDDD), normalizedElevation).lerp(colorDeviance, 0.10);
            else if (tile.biome === "mountain") terrainColor = new THREE.Color(0x555544).lerp(new THREE.Color(0x444433), normalizedElevation).lerp(colorDeviance, 0.05);
            else terrainColor = new THREE.Color(0xFF00FF);
        }
        else
        {
            var normalizedElevation = Math.min((tile.elevation - 0.8) / 0.5, 1);
            if (tile.biome === "mountain") terrainColor = new THREE.Color(0x444433).lerp(new THREE.Color(0x333322), normalizedElevation).lerp(colorDeviance, 0.05);
            else if (tile.biome === "snowyMountain") terrainColor = new THREE.Color(0xDDDDDD).lerp(new THREE.Color(0xFFFFFF), normalizedElevation).lerp(colorDeviance, 0.10);
            else terrainColor = new THREE.Color(0xFF00FF);
        }
        
        var plateColor = tile.plate.color.clone();

        var elevationColor;
        if (tile.elevation <= 0) elevationColor = new THREE.Color(0x224488).lerp(new THREE.Color(0xAADDFF), Math.max(0, Math.min((tile.elevation + 3/4) / (3/4), 1)));
        else if (tile.elevation < 0.75) elevationColor = new THREE.Color(0x997755).lerp(new THREE.Color(0x553311), Math.max(0, Math.min((tile.elevation) / (3/4), 1)));
        else elevationColor = new THREE.Color(0x553311).lerp(new THREE.Color(0x222222), Math.max(0, Math.min((tile.elevation - 3/4) / (1/2), 1)));

        var temperatureColor;
        if (tile.temperature <= 0) temperatureColor = new THREE.Color(0x0000FF).lerp(new THREE.Color(0xBBDDFF), Math.max(0, Math.min((tile.temperature + 2/3) / (2/3), 1)));
        else temperatureColor = new THREE.Color(0xFFFF00).lerp(new THREE.Color(0xFF0000), Math.max(0, Math.min((tile.temperature) / (3/3), 1)));

        var moistureColor = new THREE.Color(0xFFCC00).lerp(new THREE.Color(0x0066FF), Math.max(0, Math.min(tile.moisture, 1)));
        
        var baseIndex = planetGeometry.vertices.length;
        planetGeometry.vertices.push(tile.averagePosition);
        for (var j = 0; j < tile.corners.length; ++j)
        {
            var cornerPosition = tile.corners[j].position;
            planetGeometry.vertices.push(cornerPosition);
            planetGeometry.vertices.push(tile.averagePosition.clone().sub(cornerPosition).multiplyScalar(0.1).add(cornerPosition));

            var i0 = j * 2;
            var i1 = ((j + 1) % tile.corners.length) * 2;
            buildTileWedge(planetGeometry.faces, baseIndex, i0, i1, tile.normal);
            buildTileWedgeColors(terrainColors, terrainColor, terrainColor.clone().multiplyScalar(1-tileborder));
            buildTileWedgeColors(plateColors, plateColor, plateColor.clone().multiplyScalar(1-tileborder));
            buildTileWedgeColors(elevationColors, elevationColor, elevationColor.clone().multiplyScalar(1-tileborder));
            buildTileWedgeColors(temperatureColors, temperatureColor, temperatureColor.clone().multiplyScalar(1-tileborder));
            buildTileWedgeColors(moistureColors, moistureColor, moistureColor.clone().multiplyScalar(1-tileborder));
            for (var k = planetGeometry.faces.length - 3; k < planetGeometry.faces.length; ++k) planetGeometry.faces[k].vertexColors = terrainColors[k];
        }
    };
    
    planetGeometry.dynamic = true;
    planetGeometry.boundingSphere = new THREE.Sphere(new Vector3(0, 0, 0), 1000);
    var planetMaterial = new THREE.MeshLambertMaterial({ color: new THREE.Color(0x000000), ambient: new THREE.Color(0xFFFFFF), vertexColors: THREE.VertexColors, });
    var planetRenderObject = new THREE.Mesh(planetGeometry, planetMaterial);
    
    return{
        geometry: planetGeometry,
        terrainColors: terrainColors,
        plateColors: plateColors,
        elevationColors: elevationColors,
        temperatureColors: temperatureColors,
        moistureColors: moistureColors,
        material: planetMaterial,
        renderObject: planetRenderObject,
    };
}


function buildPlateBoundariesRenderObject(borders){
    var geometry = new THREE.Geometry();

    var i = 0;
    for(var i = 0;i<borders.length;i++)
    {
        if (i >= borders.length) return;
        
        var border = borders[i];
        if (border.betweenPlates)
        {
            var normal = border.midpoint.clone().normalize();
            var offset = normal.clone().multiplyScalar(1);

            var borderPoint0 = border.corners[0].position;
            var borderPoint1 = border.corners[1].position;
            var tilePoint0 = border.tiles[0].averagePosition;
            var tilePoint1 = border.tiles[1].averagePosition;
            
            var baseIndex = geometry.vertices.length;
            geometry.vertices.push(borderPoint0.clone().add(offset));
            geometry.vertices.push(borderPoint1.clone().add(offset));
            geometry.vertices.push(tilePoint0.clone().sub(borderPoint0).multiplyScalar(0.2).add(borderPoint0).add(offset));
            geometry.vertices.push(tilePoint0.clone().sub(borderPoint1).multiplyScalar(0.2).add(borderPoint1).add(offset));
            geometry.vertices.push(tilePoint1.clone().sub(borderPoint0).multiplyScalar(0.2).add(borderPoint0).add(offset));
            geometry.vertices.push(tilePoint1.clone().sub(borderPoint1).multiplyScalar(0.2).add(borderPoint1).add(offset));
            
            var pressure = Math.max(-1, Math.min((border.corners[0].pressure + border.corners[1].pressure) / 2, 1));
            var shear = Math.max(0, Math.min((border.corners[0].shear + border.corners[1].shear) / 2, 1));
            var innerColor = (pressure <= 0) ? new THREE.Color(1 + pressure, 1, 0) : new THREE.Color(1, 1 - pressure, 0);
            var outerColor = new THREE.Color(0, shear / 2, shear);
            
            geometry.faces.push(new THREE.Face3(baseIndex + 0, baseIndex + 1, baseIndex + 2, normal, [ innerColor, innerColor, outerColor ] ));
            geometry.faces.push(new THREE.Face3(baseIndex + 1, baseIndex + 3, baseIndex + 2, normal, [ innerColor, outerColor, outerColor ] ));
            geometry.faces.push(new THREE.Face3(baseIndex + 1, baseIndex + 0, baseIndex + 5, normal, [ innerColor, innerColor, outerColor ] ));
            geometry.faces.push(new THREE.Face3(baseIndex + 0, baseIndex + 4, baseIndex + 5, normal, [ innerColor, outerColor, outerColor ] ));
        }
        
    };
    
    geometry.boundingSphere = new THREE.Sphere(new Vector3(0, 0, 0), 1010);
    var material = new THREE.MeshBasicMaterial({ vertexColors: THREE.VertexColors, });
    var renderObject = new THREE.Mesh(geometry, material);
    
    return{
        geometry: geometry,
        material: material,
        renderObject: renderObject,
    };
}

function buildPlateMovementsRenderObject(tiles){
    var geometry = new THREE.Geometry();

    var i = 0;
    while(i<tiles.length)
    {
        if (i >= tiles.length) return;
        
        var tile = tiles[i];
        var plate = tile.plate;
        var movement = plate.calculateMovement(tile.position);
        var plateMovementColor = new THREE.Color(1 - plate.color.r, 1 - plate.color.g, 1 - plate.color.b);

        buildArrow(geometry, tile.position.clone().multiplyScalar(1.002), movement.clone().multiplyScalar(0.5), tile.position.clone().normalize(), Math.min(movement.length(), 4), plateMovementColor);
        
        tile.plateMovement = movement;
        i++
    };

    geometry.boundingSphere = new THREE.Sphere(new Vector3(0, 0, 0), 1010);
    var material = new THREE.MeshBasicMaterial({ vertexColors: THREE.VertexColors, });
    var renderObject = new THREE.Mesh(geometry, material);
    
   return{
        geometry: geometry,
        material: material,
        renderObject: renderObject,
    };
}

function buildAirCurrentsRenderObject(corners){
    var geometry = new THREE.Geometry();

    var i = 0;
    for(var i = 0;i<corners.length;i++){
        if (i >= corners.length) return;
        
        var corner = corners[i];
        buildArrow(geometry, corner.position.clone().multiplyScalar(1.002), corner.airCurrent.clone().multiplyScalar(0.5), corner.position.clone().normalize(), Math.min(corner.airCurrent.length(), 4));
    };

    geometry.boundingSphere = new THREE.Sphere(new Vector3(0, 0, 0), 1010);
    var material = new THREE.MeshBasicMaterial({ color: new THREE.Color(0xFFFFFF), });
    var renderObject = new THREE.Mesh(geometry, material);
    
    return{
        geometry: geometry,
        material: material,
        renderObject: renderObject,
    };
}

function buildArrow(geometry, position, direction, normal, baseWidth, color){
    if (direction.lengthSq() === 0) return;
    var sideOffset = direction.clone().cross(normal).setLength(baseWidth / 2);
    var baseIndex = geometry.vertices.length;
    geometry.vertices.push(position.clone().add(sideOffset), position.clone().add(direction), position.clone().sub(sideOffset));
    geometry.faces.push(new THREE.Face3(baseIndex, baseIndex + 2, baseIndex + 1, normal, [ color, color, color ]));
}

function setSurfaceRenderMode(planet, mode){

    surfaceRenderMode = mode;

    var colors;
    if (mode === "terrain") colors = planet.renderData.surface.terrainColors;
    else if (mode === "plates") colors = planet.renderData.surface.plateColors;
    else if (mode === "elevation") colors = planet.renderData.surface.elevationColors;
    else if (mode === "temperature") colors = planet.renderData.surface.temperatureColors;
    else if (mode === "moisture") colors = planet.renderData.surface.moistureColors;
    else colors = planet.renderData.surface.terrainColors;

    var faces = planet.renderData.surface.geometry.faces;
    for (var i = 0; i < faces.length; ++i)
    faces[i].vertexColors = colors[i];

    planet.renderData.surface.geometry.colorsNeedUpdate = true;
}

function showHidePlateBoundaries(show){
    if (typeof(show) === "boolean") renderPlateBoundaries = show;
    else renderPlateBoundaries = !renderPlateBoundaries;

    if (!planet) return;
    
    if (renderPlateBoundaries) planet.renderData.surface.renderObject.add(planet.renderData.plateBoundaries.renderObject);
    else planet.renderData.surface.renderObject.remove(planet.renderData.plateBoundaries.renderObject);
}

function showHidePlateMovements(show){
    if (typeof(show) === "boolean") renderPlateMovements = show;
    else renderPlateMovements = !renderPlateMovements;

    if (!planet) return;
    
    if (renderPlateMovements) planet.renderData.surface.renderObject.add(planet.renderData.plateMovements.renderObject);
    else planet.renderData.surface.renderObject.remove(planet.renderData.plateMovements.renderObject);
}

function showHideAirCurrents(show){
    if (typeof(show) === "boolean") renderAirCurrents = show;
    else renderAirCurrents = !renderAirCurrents;

    if (!planet) return;
    
    if (renderAirCurrents) planet.renderData.surface.renderObject.add(planet.renderData.airCurrents.renderObject);
    else planet.renderData.surface.renderObject.remove(planet.renderData.airCurrents.renderObject);
}



//CONTROLS
document.onkeydown = checkKey;
function checkKey(e) {

    e = e || window.event;

    // E
    if (e.keyCode == '69') {
        setSurfaceRenderMode(planet,"terrain")
    }
    //R
    else if ((e.keyCode == '82')){
        setSurfaceRenderMode(planet,"elevation")
    }

    //T
    else if ((e.keyCode == '84')){
        setSurfaceRenderMode(planet,"temperature")
    }

    //Y
    else if ((e.keyCode == '89')){
        setSurfaceRenderMode(planet,"moisture")
    }

    //U
    else if ((e.keyCode == '85')){
        setSurfaceRenderMode(planet,"plates")
    }  

    //D
    else if ((e.keyCode == '68')){
        showPlateBoundaries = !showPlateBoundaries
        showHidePlateBoundaries(showPlateBoundaries)
    } 

    //F
    else if ((e.keyCode == '70')){
        showPlateMovements = !showPlateMovements
        showHidePlateMovements(showPlateMovements)
    } 

    //G
    else if ((e.keyCode == '71')){
        showAirCurrents = !showAirCurrents
        showHideAirCurrents(showAirCurrents)
    } 

    //H
    else if ((e.keyCode == '72')){
        var selectedObject = scene.getObjectByName("planet");
        scene.remove( selectedObject );
        console.log("Extrude Planet")
        extrude(planet, 30)
        console.log("Rebuild surface")
        planet.renderData.surface = buildSurfaceRenderObject(planet.tiledTopology)


        drawPlanet(planet)
    } 

    //J
    else if ((e.keyCode == '74')){
        rotate = !rotate
    } 
    //G
    else if ((e.keyCode == '00')){
    }
    else if ((e.keyCode == '00')){
        emptyScene()
        var elevationBorderQueueSorter = function(left, right) { return left.distanceToPlateBoundary - right.distanceToPlateBoundary; };
        tiles = planet.tiledTopology.tiles;

        
        var plates = planet.plates
        plates.forEach(function(p){
            var diff = 0
            diff = p.elevation
            p.elevation = getAverageElevation(p)
            p.driftRate+=p.driftRate
            diff -=p.elevation
            p.tiles.forEach(function(t){
                t.elevation-=diff
            })
        })
        identifyBoundaryBorders(planet.tiledTopology.borders)
        var corners = collectBoundaryCorners(planet.tiledTopology.corners)
        var indexes = calculatePlateBoundaryStress(corners)


        blurPlateBoundaryStress(corners,3,0.4)

        var elevationQueue = populateElevationBorderQueue(corners, indexes)
        processElevationBorderQueue(elevationQueue, elevationBorderQueueSorter);


        calculateTileAverageElevations(planet.tiledTopology.tiles);

        setColors(planet)
        planet.geometry = generatePlanetTiledGeometry(planet);
        drawPlanet(planet);
        addLight();
        //drawBoundaryBorders(planet.plates)
        
        console.log("done")
    }
}




//HELPERS
function getEdgeOppositeFaceIndex(edge, faceIndex){
    if (edge.f[0] === faceIndex) return edge.f[1];
    if (edge.f[1] === faceIndex) return edge.f[0];
    throw "Given face is not part of given edge.";
}

function getFaceOppositeNodeIndex(face, edge){
    if (face.n[0] !== edge.n[0] && face.n[0] !== edge.n[1]) return 0;
    if (face.n[1] !== edge.n[0] && face.n[1] !== edge.n[1]) return 1;
    if (face.n[2] !== edge.n[0] && face.n[2] !== edge.n[1]) return 2;
    throw "Cannot find node of given face that is not also a node of given edge.";
}

function findNextFaceIndex(mesh, nodeIndex, faceIndex){
    var node = mesh.nodes[nodeIndex];
    var face = mesh.faces[faceIndex];
    var nodeFaceIndex = face.n.indexOf(nodeIndex);
    var edge = mesh.edges[face.e[(nodeFaceIndex + 2) % 3]];
    return getEdgeOppositeFaceIndex(edge, faceIndex);
}

function slerp(p0, p1, t){
    var omega = Math.acos(p0.dot(p1));
    return p0.clone().multiplyScalar(Math.sin((1 - t) * omega)).add(p1.clone().multiplyScalar(Math.sin(t * omega))).divideScalar(Math.sin(omega));
}

function calculateTriangleArea(pa, pb, pc){
    var vab = new THREE.Vector3().subVectors(pb, pa);
    var vac = new THREE.Vector3().subVectors(pc, pa);
    var faceNormal = new THREE.Vector3().crossVectors(vab, vac);
    var vabNormal = new THREE.Vector3().crossVectors(faceNormal, vab).normalize();
    var plane = new THREE.Plane().setFromNormalAndCoplanarPoint(vabNormal, pa);
    var height = plane.distanceToPoint(pc);
    var width = vab.length();
    var area = width * height * 0.5;
    return area;
}

function calculateFaceCentroid(pa, pb, pc){       
    var vabHalf = pb.clone().sub(pa).divideScalar(2);
    var pabHalf = pa.clone().add(vabHalf);
    var centroid = pc.clone().sub(pabHalf).multiplyScalar(1/3).add(pabHalf);
    return centroid;
}

function buildTileWedge(f, b, s, t, n){
    f.push(new THREE.Face3(b + s + 2, b + t + 2, b, n));
    f.push(new THREE.Face3(b + s + 1, b + t + 1, b + t + 2, n));
    f.push(new THREE.Face3(b + s + 1, b + t + 2, b + s + 2, n));
}

function buildTileWedgeColors(f, c, bc){
    f.push([ c, c, c ]);
    f.push([ bc, bc, c ]);
    f.push([ bc, c, c ]);
}

function extrude(planet, level){
    for(var i = 0;i<level;i++){
        planet.tiledTopology.corners.forEach(function(c){
            var normal = c.position.clone().normalize();
            var offset = normal.clone().multiplyScalar(Math.max(c.elevation,0));
            c.position.add(offset)

        })
        planet.tiledTopology.tiles.forEach(function(t){
            var normal = t.averagePosition.clone().normalize();
            var offset = normal.clone().multiplyScalar(Math.max(t.elevation,0));
            t.averagePosition.add(offset);
        })
    }
}

function getAverageElevation(plate){
    var e = 0.0
    var i = 0.0
    plate.tiles.forEach(function(t){
        e+=t.elevation
        i++
    })

    return e/i
}



//STRUCT
function Corner(id, position, cornerCount, borderCount, tileCount){
    this.id = id;
    this.position = position;
    this.corners = new Array(cornerCount);
    this.borders = new Array(borderCount);
    this.tiles = new Array(tileCount);
}

Corner.prototype.vectorTo = function Corner_vectorTo(corner){
    return corner.position.clone().sub(this.position);
};

Corner.prototype.toString = function Corner_toString(){
    return "Corner " + this.id.toFixed(0) + " < " + this.position.x.toFixed(0) + ", " + this.position.y.toFixed(0) + ", " + this.position.z.toFixed(0) + " >";
};

function Border(id, cornerCount, borderCount, tileCount){
    this.id = id;
    this.corners = new Array(cornerCount);
    this.borders = new Array(borderCount);
    this.tiles = new Array(tileCount);
}

Border.prototype.oppositeCorner = function Border_oppositeCorner(corner){
    return (this.corners[0] === corner) ? this.corners[1] : this.corners[0];
};

Border.prototype.oppositeTile = function Border_oppositeTile(tile){
    return (this.tiles[0] === tile) ? this.tiles[1] : this.tiles[0];
};

Border.prototype.length = function Border_length(){
    return this.corners[0].position.distanceTo(this.corners[1].position);
};

Border.prototype.isLandBoundary = function Border_isLandBoundary(){
    return (this.tiles[0].elevation > 0) !== (this.tiles[1].elevation > 0);
};

Border.prototype.toString = function Border_toString(){
    return "Border " + this.id.toFixed(0);
};

function Tile(id, position, cornerCount, borderCount, tileCount){
    this.id = id;
    this.position = position;
    this.corners = new Array(cornerCount);
    this.borders = new Array(borderCount);
    this.tiles = new Array(tileCount);
}

function Plate(color, driftAxis, driftRate, spinRate, elevation, oceanic, root){
    this.color = color;
    this.driftAxis = driftAxis;
    this.driftRate = driftRate;
    this.spinRate = spinRate;
    this.elevation = elevation;
    this.oceanic = oceanic;
    this.root = root;
    this.tiles = [];
    this.boundaryCorners = [];
    this.boundaryBorders = [];
}

Plate.prototype.calculateMovement = function Plate_calculateMovement(position){
    var movement = this.driftAxis.clone().cross(position).setLength(this.driftRate * position.clone().projectOnVector(this.driftAxis).distanceTo(position));
    movement.add(this.root.position.clone().cross(position).setLength(this.spinRate * position.clone().projectOnVector(this.root.position).distanceTo(position)));
    return movement.multiplyScalar(1);
};



//RANDOM
function randomFromInterval(min,max){
    return Math.random()*(max-min)+min;
}

function randomIntFromInterval(min,max){
    return Math.floor(Math.random()*(max-min+1)+min);
}

function randomUnitVector(random){
    var theta = Math.floor(randomFromInterval(0, Math.PI * 2));
    var phi = Math.acos(Math.floor(randomFromInterval(-1, 1)));
    var sinPhi = Math.sin(phi);
    return new Vector3(
        Math.cos(theta) * sinPhi,
        Math.sin(theta) * sinPhi,
        Math.cos(phi));
}

randomReal = function (min, max){
    return randomUnit() * (max - min) + min;
};

randomRealInclusive = function(min, max){
    return randomUnitInclusive() * (max - min) + min;
};

randomUnit = function (){
    return this.randomNext() / 0x80000000;
};
randomUnitInclusive = function (){
    return randomNext() / 0x7FFFFFFF;
};

randomNext = function(){
    return Math.random();
};

randomInteger = function (min, max){
    return this.randomIntegerExclusive(min, max + 1);
};

randomIntegerExclusive = function (min, max){
    min = Math.floor(min);
    max = Math.floor(max);
    return Math.floor(this.randomUnit() * (max - min)) + min;
};

