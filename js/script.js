    
var scene;
var camera;
var renderer;
var projector;
var controls;
var light;
// array of functions for the rendering loop
var onRenderFcts = [];
var Vector3 = THREE.Vector3
var orderLevel = 2

var planet = {};
var planetMesh;
var wireframe = false;


window.onload = function(){
    initializeScene();
    generatePlanet();
    renderScene();      
}

function initializeScene(){
    // init renderer
    renderer    = new THREE.WebGLRenderer({
        antialias   : true
    });
    renderer.setClearColor(new THREE.Color('lightgrey'), 1)
    renderer.setFaceCulling(THREE.CullFaceFront, THREE.FrontFaceDirectionCW);
    

    renderer.setSize( window.innerWidth, window.innerHeight );
    document.body.appendChild( renderer.domElement );

    projector = new THREE.Projector();



    // init scene and camera
    scene = new THREE.Scene();
    
    camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.01, 1000);
    camera.position.z = 5;

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




function generatePlanet(){

    //generate subdivided icosahedron
    planet.topology = generateSubdividedIcosahedron(orderLevel);

    // distort mesh
    rotateEdge(planet.topology,89)

    // calculate centroids 
    centroids(planet.topology)

    //reorder triangle nodes
    reorderTriangleNodes(planet.topology)


    //compute dual polyhedron
    planet.tiledTopology = generatePlanetTiledTopology(planet.topology);

    //
    setRandomColors()

    //compute dual poly geometry
    planet.geometry = generatePlanetTiledGeometry();
    
    //create mesh
    planetMesh = generatePlanetMesh();

    scene.add(planetMesh)

    console.log("Planet: ")      
    console.log(planet)  
}

//Topology 
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


function generatePlanetTiledTopology(topology){
    var corners = new Array(topology.faces.length);
    var borders = new Array(topology.edges.length);
    var tiles = new Array(topology.nodes.length);
    

    for (var i = 0; i < topology.faces.length; ++i){
        var face = topology.faces[i];
        corners[i] = new Corner(i, face.centroid.clone().multiplyScalar(1), face.e.length, face.e.length, face.n.length);
    }

    for (var i = 0; i < topology.edges.length; ++i){
        var edge = topology.edges[i];
        borders[i] = new Border(i, 2, 4, 2); //edge.f.length, topology.faces[edge.f[0]].e.length + topology.faces[edge.f[1]].e.length - 2, edge.n.length
    }

    for (var i = 0; i < topology.nodes.length; ++i){
        var node = topology.nodes[i];
        tiles[i] = new Tile(i, node.p.clone().multiplyScalar(1), node.f.length, node.e.length, node.e.length);
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
            if(tile.borders[j] == undefined)
                console.log(i)
            //area += calculateTriangleArea(tile.position, tile.borders[j].corners[0].position, tile.borders[j].corners[1].position);
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

//
function rotateEdge(mesh, edgeIndex){

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

    var predicate = function(oldNode0, oldNode1, newNode0, newNode1){
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

    if (!predicate(oldNode0, oldNode1, newNode0, newNode1)) return false;
    
    oldNode0.e.splice(oldNode1.e.indexOf(edgeIndex), 1);

    oldNode1.e.splice(oldNode0.e.indexOf(edgeIndex), 1);
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
    newNode0.f.push(edge.f[0]);
    newNode1.f.push(edge.f[1]);
    
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



//Geometry
function generatePlanetGeometry(){
    
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

function generatePlanetTiledGeometry(){
    
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


//Mesh
function generatePlanetMesh(){
    var material    = new THREE.MeshBasicMaterial({
        color:0xff0000,
        wireframe: true,
        side:THREE.DoubleSide 
    }); 
    var material = new THREE.MeshNormalMaterial();
    var material = new THREE.MeshLambertMaterial({ color: new THREE.Color(0x000000), ambient: new THREE.Color(0xFFFFFF), vertexColors: THREE.VertexColors, });

    var mesh    = new THREE.Mesh(planet.geometry ,material);

    return mesh;
}




//color
function modifyTileColor(i, color){
    planet.tiledTopology.tiles[i].color = color
}

function setRandomColors(){
    var colorDeviance = new THREE.Color(0.2,0.2,0.2);
    for(var i = 0; i < planet.tiledTopology.tiles.length;i++){
        planet.tiledTopology.tiles[i].elevation = -Math.random()
    }

    for(var i = 0; i < planet.tiledTopology.tiles.length;i++){
        var tile = planet.tiledTopology.tiles[i]
        modifyTileColor(i, new THREE.Color(0x0066FF).lerp(new THREE.Color(0x0044BB), Math.min(-tile.elevation, 1)).lerp(colorDeviance, 0.10))
        modifyTileColor(i, new THREE.Color('#'+Math.floor(Math.random()*16777215).toString(16)))
    }
}

function buildTileWedge(f, b, s, t, n){
    f.push(new THREE.Face3(b + s + 2, b + t + 2, b, n));
    f.push(new THREE.Face3(b + s + 1, b + t + 1, b + t + 2, n));
    f.push(new THREE.Face3(b + s + 1, b + t + 2, b + s + 2, n));
}







document.onkeydown = checkKey;

function checkKey(e) {

    e = e || window.event;


    if (e.keyCode == '81') {        
        for( var i = scene.children.length - 1; i >= 0; i--) {obj = scene.children[i];scene.remove(obj);}
        
        planet.geometry = generatePlanetGeometry();
        scene.add(generatePlanetMesh())
        addLight();
    }
    else if ((e.keyCode == '65')){
        for( var i = scene.children.length - 1; i >= 0; i--) {obj = scene.children[i];scene.remove(obj);}
        //rotateEdge(planet.topology,2)
        planet.tiledTopology = generatePlanetTiledTopology(planet.topology); 
        setRandomColors()
        planet.geometry = generatePlanetTiledGeometry();
        scene.add(generatePlanetMesh())
        addLight();
    }
    else if ((e.keyCode == '87')){
        rotateEdge(planet.topology,89)

    }
    else if ((e.keyCode == '90')){
        if(wireframe == false){
            var material    = new THREE.MeshBasicMaterial({
                color:0xff0000,
                wireframe: true,
                side:THREE.DoubleSide 
            }); 
            var mesh    = new THREE.Mesh(planet.geometry   ,material);
            scene.add(mesh)
        }
            
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

function calculateCentroids(geometry){
    for (var i = 0; i < geometry.faces.length; ++i)
    {
        var face = geometry.faces[i];
        var p0 = geometry.vertices[face.a]
        var p1 = geometry.vertices[face.b]
        var p2 = geometry.vertices[face.c]
        face.centroid = calculateFaceCentroid(p0, p1, p2);
    }
}













//FIRST TEST
function generateIcoGeometry(){
    var phi = (1.0 + Math.sqrt(5.0)) / 2.0;

    var geometry = new THREE.Geometry();

    geometry.vertices.push(new THREE.Vector3( -1,  phi, 0 ).normalize());
    geometry.vertices.push(new THREE.Vector3( 1,  phi, 0 ).normalize());
    geometry.vertices.push(new THREE.Vector3( -1,  -phi, 0 ).normalize());
    geometry.vertices.push(new THREE.Vector3( 1,  -phi, 0 ).normalize());

    geometry.vertices.push(new THREE.Vector3( 0,  -1, phi ).normalize());
    geometry.vertices.push(new THREE.Vector3( 0,  1, phi ).normalize());
    geometry.vertices.push(new THREE.Vector3( 0,  -1, -phi ).normalize());
    geometry.vertices.push(new THREE.Vector3( 0,  1, -phi ).normalize());

    geometry.vertices.push(new THREE.Vector3( phi,  0, -1 ).normalize());
    geometry.vertices.push(new THREE.Vector3( phi,  0, 1 ).normalize());
    geometry.vertices.push(new THREE.Vector3( -phi,  0, -1 ).normalize());
    geometry.vertices.push(new THREE.Vector3( -phi,  0, 1 ).normalize());


    
    geometry.faces.push(new THREE.Face3( 0,  11, 5 ));
    geometry.faces.push(new THREE.Face3( 0,  5, 1 ));
    geometry.faces.push(new THREE.Face3( 0,  1, 7 ));
    geometry.faces.push(new THREE.Face3( 0,  7, 10 ));
    geometry.faces.push(new THREE.Face3( 0,  10, 11 ));

    geometry.faces.push(new THREE.Face3(1, 5, 9));
    geometry.faces.push(new THREE.Face3(5, 11, 4));
    geometry.faces.push(new THREE.Face3(11, 10, 2));
    geometry.faces.push(new THREE.Face3(10, 7, 6));
    geometry.faces.push(new THREE.Face3(7, 1, 8));

    geometry.faces.push(new THREE.Face3(3, 9, 4));
    geometry.faces.push(new THREE.Face3(3, 4, 2));
    geometry.faces.push(new THREE.Face3(3, 2, 6));
    geometry.faces.push(new THREE.Face3(3, 6, 8));
    geometry.faces.push(new THREE.Face3(3, 8, 9));

    geometry.faces.push(new THREE.Face3(4, 9, 5));
    geometry.faces.push(new THREE.Face3(2, 4, 11));
    geometry.faces.push(new THREE.Face3(6, 2, 10));
    geometry.faces.push(new THREE.Face3(8, 6, 7));
    geometry.faces.push(new THREE.Face3(9, 8, 1));

    geometry.computeBoundingSphere();
    geometry.computeVertexNormals();
    geometry.computeFaceNormals();

    return geometry;
}

function generateSubdividedIcoGeometry(order){
    var geometry = generateIcoGeometry();

    for (var i = 0; i < order; i++)
    {
        var geometry2 = new THREE.Geometry();
        geometry2.vertices = geometry.vertices.slice(0);
        var faces2 = [];

        var faceNb = 0
        geometry.faces.forEach(function(f){ 

            var a = getMiddlePoint(geometry2,geometry2.vertices[f.a],geometry2.vertices[f.b]);
            var b = getMiddlePoint(geometry2,geometry2.vertices[f.b],geometry2.vertices[f.c]);
            var c = getMiddlePoint(geometry2,geometry2.vertices[f.c],geometry2.vertices[f.a]); 

            geometry2.faces.push(new THREE.Face3(indexOfVertice(geometry2.vertices,a), indexOfVertice(geometry2.vertices,b), indexOfVertice(geometry2.vertices,c)));
            geometry2.faces.push(new THREE.Face3(indexOfVertice(geometry2.vertices,a), indexOfVertice(geometry2.vertices,c), f.a));
            geometry2.faces.push(new THREE.Face3(indexOfVertice(geometry2.vertices,b), indexOfVertice(geometry2.vertices,a), f.b));
            geometry2.faces.push(new THREE.Face3(indexOfVertice(geometry2.vertices,c), indexOfVertice(geometry2.vertices,b), f.c));
        })

        geometry = geometry2
    }

    geometry.computeBoundingSphere();
    geometry.computeVertexNormals();
    geometry.computeFaceNormals();

    //console.log("Geometry : ")
    //console.log(geometry)
    return geometry;
}

function generateTiledGeometry(order){
    var geometry = generateSubdividedIcoGeometry(order);
    calculateCentroids(geometry)
    var geometry2 = new THREE.Geometry();

    for (var i = 0;i<geometry.vertices.length;i++){
        tryAddVertice(geometry2,geometry.vertices[i])
    }
    for (var i = 0;i<geometry.faces.length;i++){
        tryAddVertice(geometry2,geometry.faces[i].centroid)
    }

    

    

    //console.log(geometry2)
}

function generateMesh(order){
    var geometry    = generateSubdividedIcoGeometry(order);
    var material    = new THREE.MeshBasicMaterial({
        color:0xff0000,
        wireframe: true
    }); 
    material = new THREE.MeshNormalMaterial();

    var mesh  = new THREE.Mesh(geometry,material);
    //console.log("Mesh :" )
    //console.log(mesh)
    //console.log("---")
    return mesh;
}

function tryAddVertice(g, vertice){

    var test = 0;
    for (var i = 0; i< g.vertices.length;i++){
        if (g.vertices[i].x == vertice.x && g.vertices[i].y == vertice.y && g.vertices[i].z == vertice.z)
            test+=1;
    }
    if (test == 0)
        g.vertices.push(vertice)

    return vertice;
}

Array.prototype.contains = function(obj) {
    var i = this.length;
    while (i--) {
        if (this[i] === obj) {
            return true;
        }
    }
    return false;
}

function calculateFaceCentroid(a, b, c){       
    var centroid = new THREE.Vector3((a.x+b.x+c.x)/3,(a.y+b.y+c.y)/3,(a.z+b.z+c.z)/3);
    return centroid;
}

function indexOfVertice(vertices, v){
    for (var i = 0; i< vertices.length;i++){
        if (vertices[i].x == v.x && vertices[i].y == v.y && vertices[i].z == v.z)
            return i;
    }
    return -1
}

function getMiddlePoint(g,p1, p2){
    var middle = new THREE.Vector3(
        (p1.x + p2.x) / 2.0, 
        (p1.y + p2.y) / 2.0, 
        (p1.z + p2.z) / 2.0).normalize();



    var test = 0;
    for (var i = 0; i< g.vertices.length;i++){
        if (g.vertices[i].x == middle.x && g.vertices[i].y == middle.y && g.vertices[i].z == middle.z)
            test+=1;
    }
    if (test == 0)
        g.vertices.push(middle)

    return middle;
}



