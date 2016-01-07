    
var scene;
var camera;
var renderer;
var projector;
var controls;
var light;
// array of functions for the rendering loop
var onRenderFcts = [];
var Vector3 = THREE.Vector3

var planet = {};
var planetMesh;
var wireframe = false;
var mode = "tile"
var renderMode = "elevation"
var noise


window.onload = function(){
    initializeScene();

    planet = generatePlanet(20,0.5);
    drawPlanet(planet)
    //drawEdges(planet)
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
    })
}


function generatePlanet(orderLevel,topologyDistortionRate){
    var planet = {};
    planet.material = [new THREE.MeshLambertMaterial({ color: new THREE.Color(0x000000), ambient: new THREE.Color(0xFFFFFF), vertexColors: THREE.VertexColors, })];
    //generate subdivided icosahedron
    planet.topology = generateSubdividedIcosahedron(orderLevel);


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
        console.log("relaxed")
    }
    
    // calculate centroids 
    centroids(planet.topology)

    //reorder triangle nodes
    reorderTriangleNodes(planet.topology)


    //compute dual polyhedron
    planet.tiledTopology = generatePlanetTiledTopology(planet.topology);

    //
    setElevation(planet);
    setColors(planet)

    //compute dual poly geometry
    planet.geometry = generatePlanetTiledGeometry(planet);

    return planet;
}


function drawPlanet(planet){
    //create mesh
    planet.mesh = generatePlanetMesh(planet);
    scene.add(planet.mesh)
}

function drawEdges(planet){
    var material = new THREE.LineBasicMaterial({
        color: 0x0000ff
    });
    var geometry = new THREE.Geometry();
    var line = new THREE.Line(geometry, material);

    planet.tiledTopology.borders.forEach(function(b){
        var a = b.corners[0].position;
        var b = b.corners[1].position;
        geometry.vertices.push(a);
        geometry.vertices.push(b); 
        var line = new THREE.Line(geometry, material);
        scene.add(line);
        geometry = new THREE.Geometry();
    })  
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

    console.log("Changed edge")
    //console.log(planet)
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

    planet.tiledTopology = { corners: corners, borders: borders, tiles: tiles }
    
   return  { corners: corners, borders: borders, tiles: tiles };
}


//GEOMETRY
function generatePlanetGeometry(planet){
    
    var mesh = planet.topology  
    var geometry = new THREE.Geometry();

    for (var i = 0; i< mesh.nodes.length; i++){
        geometry.vertices.push(mesh.nodes[i].p);
    }

    for (var i = 0; i< mesh.faces.length; i++){
        geometry.faces.push(new THREE.Face3( mesh.faces[i].n[2],   mesh.faces[i].n[0],mesh.faces[i].n[1] ));
    }

    geometry.computeBoundingSphere();
    geometry.computeVertexNormals();
    geometry.computeFaceNormals();

    return geometry;
}

function generatePlanetTiledGeometry(planet){
    
    var geometry = new THREE.Geometry();
    var tiles = planet.tiledTopology.tiles

    for(var i = 0;i<tiles.length;i++){

        var tile = tiles[i];
        var baseIndex = geometry.vertices.length;
        geometry.vertices.push(tile.averagePosition);
        for (var j = 0; j < tile.corners.length; ++j)
        {
            var cornerPosition = tile.corners[j].position;
            geometry.vertices.push(cornerPosition);
            geometry.vertices.push(tile.averagePosition.clone().sub(cornerPosition).multiplyScalar(1).add(cornerPosition));

            var i0 = j * 2;
            var i1 = ((j + 1) % tile.corners.length) * 2;

            buildTileWedge(geometry.faces, baseIndex, i0, i1, tile.normal)
            for (var k = geometry.faces.length - 3; k < geometry.faces.length; ++k){
                if (tile.color!= undefined)geometry.faces[k].color = tile.color;
            } 
        }
    }

    planet.geometry = geometry;

    return geometry;
}


//MESH
function generatePlanetMesh(planet){
    planet.mesh = new THREE.SceneUtils.createMultiMaterialObject(planet.geometry ,planet.material);
    return planet.mesh;
}




//color
function modifyTileColor(i, color){
    planet.tiledTopology.tiles[i].color = color
}

function setColors(planet){
    var colorDeviance = new THREE.Color(Math.random(), Math.random(), Math.random());

    for(var i = 0; i < planet.tiledTopology.tiles.length;i++){
        var tile = planet.tiledTopology.tiles[i]
        if(renderMode == "biome"){
            if(tile.elevation<0)
            modifyTileColor(i, new THREE.Color(0x0066FF).lerp(new THREE.Color(0x0044BB), Math.min(-tile.elevation, 1)).lerp(colorDeviance, 0.10))
        else
            modifyTileColor(i, new THREE.Color(0xAA9977).lerp(new THREE.Color(0x887755), tile.elevation).lerp(colorDeviance, 0.10))
        }
        
        else if (renderMode == "elevation"){
            if (tile.elevation <= 0)
            modifyTileColor(i, new THREE.Color(0x224488).lerp(new THREE.Color(0xAADDFF), Math.max(0, Math.min((tile.elevation + 3/4) / (3/4), 1))))

        else if (tile.elevation < 0.75)
            modifyTileColor(i, new THREE.Color(0x997755).lerp(new THREE.Color(0x553311), Math.max(0, Math.min((tile.elevation) / (3/4), 1))))

        else
            modifyTileColor(i,  new THREE.Color(0x553311).lerp(new THREE.Color(0x222222), Math.max(0, Math.min((tile.elevation - 3/4) / (1/2), 1))))
        }
        else if (renderMode == "plates"){
            modifyTileColor(i, tile.plate.color)
        }

    }
}

function setElevation(planet){
    tiles = planet.tiledTopology.tiles;
    generatePlanetTerrain(planet)

    //for(var i = 0; i < planet.tiledTopology.tiles.length;i++){
    //    x = planet.tiledTopology.tiles[i].position.x;
    //   y = planet.tiledTopology.tiles[i].position.y;
    //    z = planet.tiledTopology.tiles[i].position.z;
    //    planet.tiledTopology.tiles[i].elevation = noise.simplex3(x , y , z)-0.3;
    //}
}

function generatePlanetTerrain(planet){

    var elevationBorderQueueSorter = function(left, right) { return left.distanceToPlateBoundary - right.distanceToPlateBoundary; };
    tiles = planet.tiledTopology.tiles;

    var plates = generateTectonicPlates(planet,planet.tiledTopology,50,0.7)
    identifyBoundaryBorders(planet.tiledTopology.borders)
    var corners = collectBoundaryCorners(planet.tiledTopology.corners)
    var indexes = calculatePlateBoundaryStress(corners)


    blurPlateBoundaryStress(corners,3,0.4)

    var elevationQueue = populateElevationBorderQueue(corners, indexes)
    processElevationBorderQueue(elevationQueue, elevationBorderQueueSorter);

    //planet.tiledTopology.corners.forEach(function(c){
    //    c.elevation-=2;
    //})

    calculateTileAverageElevations(planet.tiledTopology.tiles);
    for(var i = 0;i<30;i++)
        extrude(planet)

    //drawBoundaryBorders(plates)
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
    
    console.log(plates)
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

function extrude(planet){
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




document.onkeydown = checkKey;
function checkKey(e) {

    e = e || window.event;

    // E
    if (e.keyCode == '69') {
        renderMode = "elevation"
        emptyScene()
        setColors(planet)
        planet.geometry = generatePlanetTiledGeometry(planet);
        drawPlanet(planet);
        
        addLight();
        
    }
    //R
    else if ((e.keyCode == '82')){
        renderMode = "biome"
        emptyScene()
        setColors(planet)
        planet.geometry = generatePlanetTiledGeometry(planet);
        drawPlanet(planet);
        addLight();
    }

    //T
    else if ((e.keyCode == '84')){
        renderMode = "plates"
        emptyScene()
        setColors(planet)
        planet.geometry = generatePlanetTiledGeometry(planet);
        drawPlanet(planet);
        addLight();
    }
    //G
    else if ((e.keyCode == '71')){
        //drawBoundaryBorders(planet.plates)
        emptyScene()
        extrude(planet)
        setColors(planet)
        planet.geometry = generatePlanetTiledGeometry(planet);
        drawPlanet(planet);
        addLight();
    }
    //Y
    else if ((e.keyCode == '89')){
        emptyScene()
        var elevationBorderQueueSorter = function(left, right) { return left.distanceToPlateBoundary - right.distanceToPlateBoundary; };
        tiles = planet.tiledTopology.tiles;

        var plates = planet.plates
        plates.forEach(function(p){
            p.elevation = getAverageElevation(p)
            p.driftRate*=2
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

function getAverageElevation(plate){
    var e = 0.0
    var i = 0.0
    plate.tiles.forEach(function(t){
        e+=t.elevation
        i++
    })
    plate.elevation = e/i
    return e/i
}

function getMaterial(planet){
    if(wireframe){
        mode == "triangle"? material = materialTriangleWireframe : material = materialTileWireframe
    }
    else{
        mode == "triangle"? material = materialTriangle : material = materialTile
    }
    return material
}

var materialTriangle =          [new THREE.MeshNormalMaterial()];
var materialTriangleWireframe = [new THREE.MeshNormalMaterial(),
                                 new THREE.MeshBasicMaterial( { color: 0x000000, shading: THREE.FlatShading, wireframe: true } ) 
                                ];
var materialTile =              [new THREE.MeshLambertMaterial({ color: new THREE.Color(0x000000), ambient: new THREE.Color(0xFFFFFF), vertexColors: THREE.VertexColors, })];
var materialTileWireframe =     [new THREE.MeshLambertMaterial({ color: new THREE.Color(0x000000), ambient: new THREE.Color(0xFFFFFF), vertexColors: THREE.VertexColors, }),
                                 new THREE.MeshBasicMaterial( { color: 0x000000, shading: THREE.FlatShading, wireframe: true } ) 
                                ];


































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

function randomFromInterval(min,max)
{
    return Math.random()*(max-min)+min;
}

function randomUnitVector(random)
{
    var theta = Math.floor(randomFromInterval(0, Math.PI * 2));
    var phi = Math.acos(Math.floor(randomFromInterval(-1, 1)));
    var sinPhi = Math.sin(phi);
    return new Vector3(
        Math.cos(theta) * sinPhi,
        Math.sin(theta) * sinPhi,
        Math.cos(phi));
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

function Plate(color, driftAxis, driftRate, spinRate, elevation, oceanic, root)
{
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


//
function emptyScene(planet){
    for( var i = scene.children.length - 1; i >= 0; i--){
        obj = scene.children[i];
        scene.remove(obj);
    }
}







