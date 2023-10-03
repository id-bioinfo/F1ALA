const width = window.innerWidth;
const height = window.innerHeight;
const rootX = width / 2;
const rootY = 100
const color = d3.scaleOrdinal(d3.schemeSet3);

function isRoot(node) {
    return node.isRoot
};

function isAnnotated(node) {
	return node.hasOwnProperty("annotation")
}

const all_links = data.links;
let links = null;

const all_nodes = sortByNumberOfSamples(generateNodesFromJSON(data))
function generateNodesFromJSON(jsonData) {
    function addFixedPositionToRoot(node) {
        return isRoot(node) ?
            { ...node, fx: rootX, fy: rootY } :
            { ...node }
    }
    return jsonData.nodes.map(addFixedPositionToRoot);
}
function sortByNumberOfSamples(nodes) {
    return nodes.sort((a, b) => a.numberOfSamples - b.numberOfSamples)
}
let nodes = null

let simulation = null;
function generateSimulation(nodes, links) {
    return d3.forceSimulation(nodes)
        .force("link", d3.forceLink(links)
            .distance(link => linkDistanceScale(link.distance))
            .id(link => link.id)
			.strength(1.5))
        .force("charge", d3.forceManyBody().strength(-200))
        .force("collision", d3.forceCollide(26))
        .force("position", d3.forceY().y(10000).strength(0.002))
        .on("tick", draw);
}

const largestLinkDistancePx = 500;
const smallestLinkDistancePx = 180;
const linkDistanceScale = d3.scaleLinear()
    .domain([d3.min(all_links, d => d.distance), d3.max(all_links, d => d.distance)])
    .range([smallestLinkDistancePx, largestLinkDistancePx])

const totalNumberOfNodes = data.nodes.length;
const sliderValueToNumberOfNodesScale = d3.scaleLinear()
    .domain([1, 1000])
    .range([1, totalNumberOfNodes])
    .clamp(true);

const defaultNumberOfNodes = 50;
initialize()
function initialize() {
    d3.select("#numberOfNodesThreshold")
        .html(defaultNumberOfNodes);
    d3.select("#numberOfNodesSlider")
        .property("value", sliderValueToNumberOfNodesScale.invert(defaultNumberOfNodes));
    nodes = findTopNNodes(all_nodes, defaultNumberOfNodes);
    links = recontructLinks(all_nodes, all_links, defaultNumberOfNodes);
    simulation = generateSimulation(nodes, links);
}

const dpi = window.devicePixelRatio;
const canvas = d3.select("#graph")
    .append("canvas")
    .attr("width", dpi * width)
    .attr("height", dpi * height)
    .attr("style", `width: ${width}px; max-width: 100%; height: auto;`)
    .node();
const context = canvas.getContext("2d");
context.scale(dpi, dpi);

function draw() {
    context.save();
    context.clearRect(0, 0, width, height);
    context.translate(transform.x, transform.y);
    context.scale(transform.k, transform.k);
    context.globalAlpha = 0.6;
    context.strokeStyle = "#999";
    context.lineWidth = 3;
    context.beginPath();
    links.forEach(drawLink);
    context.stroke();
    context.strokeStyle = "#fff";
    context.globalAlpha = 1;
    nodes.forEach(node => {
        context.beginPath();
        drawNode(node)
        context.fillStyle = isAnnotated(node) ? color(node.annotation) : "#ccc";
        context.fill();
        context.stroke();
    });
    context.fillStyle = "#000";
    context.textAlign = "center";
    context.textBaseline = 'middle';
    context.font = "15px Arial";
    nodes.forEach(node => {
        drawNodeText(node)
    });
    context.restore();
}

const linkArrowHeadLength = 30;
function drawLink(d) {

    context.moveTo(d.source.x, d.source.y);
    context.lineTo(d.target.x, d.target.y);
}

function drawNode(node) {
    const nodeRadius = nodeRadiusScale(node.numberOfSamples)
    context.moveTo(node.x + nodeRadius, node.y);
    context.arc(node.x, node.y, nodeRadius, 0, 2 * Math.PI);
}

function drawNodeText(node) {
	const nodeLabel = constructNodeLabel(node);
    context.fillText(nodeLabel, node.x, node.y - 10);
    const numberOfSamplesLabel = node.numberOfSamples;
    context.fillText(numberOfSamplesLabel, node.x, node.y + 10);
}

function constructNodeLabel(node) {
	const annotationLabel = isAnnotated(node) ? node.annotation : ""
	const rootLabel = isRoot(node) ? "[Root]" : ""
	return `${rootLabel}${annotationLabel}`;
}

const largestNodeRadiusPx = 100;
const smallestNodeRaidusPx = 25;
const nodeRadiusScale = d3.scaleSqrt()
    .domain([d3.min(all_nodes, n => n.numberOfSamples), d3.max(all_nodes, n => n.numberOfSamples)])
    .range([smallestNodeRaidusPx, largestNodeRadiusPx])

let transform = d3.zoomIdentity;

d3.select(canvas).call(
    d3.drag()
        .subject(findSubject)
        .on("start", dragStarted)
        .on("drag", dragged)
        .on("end", dragEnded)
)

function findSubject(event) {
    const px = transform.invertX(event.x);

    const py = transform.invertY(event.y);
    const selectedNode = findNearestNode(px, py)
    if (selectedNode && !isRoot(selectedNode)) {
        selectedNode.x = transform.applyX(selectedNode.x);
        selectedNode.y = transform.applyY(selectedNode.y);
        return selectedNode
    }
}

function findNearestNode(px, py) {
    return d3.least(nodes, getNodeDistance({ px, py }));
}

const getNodeDistance = ({ px, py }) => ({ x, y, numberOfSamples }) => {
    const pointerToNodeDistance = (x - px) ** 2 + (y - py) ** 2;
    const nodeRadius = (nodeRadiusScale(numberOfSamples) + 20) ** 2;
    if (pointerToNodeDistance < nodeRadius) {
        return pointerToNodeDistance;
    }
}

function dragStarted(event) {
    const node = event.subject;
    if (!isRoot(node)) {
        if (!event.active) {
            simulation.alphaTarget(0.3).restart();
        }
        node.fx = transform.invertX(event.x);
        node.fy = transform.invertY(event.y);
    }
}

function dragged(event) {
    const node = event.subject;
    if (!isRoot(node)) {
        node.fx = transform.invertX(event.x);
        node.fy = transform.invertY(event.y);
    }
}

function dragEnded(event) {
    const node = event.subject;
    if (!isRoot(node)) {
        if (!event.active) {
            simulation.alphaTarget(0);
        }
        node.fx = null;
        node.fy = null;
    }
}

d3.select(canvas).call(
    d3.zoom()
        .scaleExtent([1 / 10, 8])
        .on("zoom", zoomed)
);

function zoomed(event) {
    transform = event.transform;
    draw();
}

d3.select(canvas).on("mousemove", mousemoved);

function mousemoved(event) {
    const px = transform.invertX(event.x);
    const py = transform.invertY(event.y);
    const hoveredNode = findNearestNode(px, py)
    if (hoveredNode) {
        showTooltip(event.x, event.y, hoveredNode)
    } else {
        hideTooltip()
    }
}

function showTooltip(x, y, node) {
    d3.select("#tooltip")
        .style("display", "block")
        .style("top", y + 15 + "px")
        .style("left", x + 15 + "px")
    setTooltipData(node)
}

function setTooltipData(node) {
    if (!isAnnotated(node)) {
        d3.select("#tooltip")
            .html(null)
            .append("p")
            .html(`Total number of samples: ${node.numberOfSamples}`);
    } else {
		let tooltipHtml = ""
		if (!isRoot(node)) {
			tooltipHtml += `Annotated Node: ${node.cladeRoot}<br/>\nDistance to root: ${node.distanceToRoot}<br/>\n`
		}
        tooltipHtml += `
        ${annotationTooltipName}: ${node.annotation}<br/>
        F1 score (annotation confidence): ${node.f1score}<br/>
        True positive samples: ${node.numberOfTruePositive}<br/>
        False positive samples: ${node.numberOfFalsePositive}<br/>
        False negative samples: ${node.numberOfFalseNegative}<br/>
        Unlabeled samples: ${node.numberOfUnlabeledSamples}
        `
        d3.select("#tooltip")
            .html(null)
            .append("p")
            .html(tooltipHtml);
    }
}

function hideTooltip() {
    d3.select("#tooltip")
        .style("display", "none");
}

d3.select("#numberOfNodesSlider")
    .on("input", sliderMoved)
    .on("mouseup", slideUped)
    .on("touchend", slideUped);

function sliderMoved(event) {
    const sliderValue = event.target.value
    const numberOfNodesThreshold = Math.ceil(sliderValueToNumberOfNodesScale(sliderValue));
    d3.select("#numberOfNodesThreshold").html(numberOfNodesThreshold);
}

function slideUped(event) {
    const sliderValue = event.target.value
    const numberOfNodesThreshold = Math.ceil(sliderValueToNumberOfNodesScale(sliderValue));
    d3.select("#numberOfNodesThreshold").html(numberOfNodesThreshold);
    nodes = findTopNNodes(all_nodes, numberOfNodesThreshold);
    links = recontructLinks(all_nodes, all_links, numberOfNodesThreshold);
    simulation.stop();
    simulation = generateSimulation(nodes, links);
    draw();
}

function findTopNNodes(all_nodes, n) {
    return all_nodes.slice(all_nodes.length - n)
}

function recontructLinks(all_nodes, all_links, numberOfNodesThreshold) {
    const deletedNodesIds = all_nodes.slice(0, all_nodes.length - numberOfNodesThreshold).map(n => n.id);
    const deletedNodesParentJSON = {};
    const output_links = [];
    for (const link of all_links) {
        if (deletedNodesIds.includes(link.target)) {
            deletedNodesParentJSON[link.target] = link.source; 
        } else {
            output_links.push(Object.assign({}, link));
        }
    }
    for (const link of output_links) {
        if (deletedNodesIds.includes(link.source)) {
            link.source = findUndeletedParent(link.source, deletedNodesParentJSON);
        }
    }
    return output_links;
}

function findUndeletedParent(nodeId, deletedNodesParentJSON) {
    if (nodeId in deletedNodesParentJSON) {
        return findUndeletedParent(deletedNodesParentJSON[nodeId], deletedNodesParentJSON);
    } else {
        return nodeId
    }
}


const saveAsPNGButton = document.getElementById("saveAsPNGButton");
saveAsPNGButton.addEventListener("click", function() {
	var image = canvas.toDataURL("image/png", 1.0)//.replace("image/png", "image/octet-stream");
	const createEl = document.createElement('a');
    createEl.href = image;

    // This is the name of our downloaded file
    createEl.download = "graph.png";
	
	//const { jsPDF } = require("jspdf");
	//var PDF = new JsPDF();
	//PDF.addImage(image, 'png', 0, 0, canvas.width, canvas.height)
	//PDF.save('graph.pdf');
    // Click the download button, causing a download, and then remove it
    createEl.click();
    createEl.remove();
}, false);