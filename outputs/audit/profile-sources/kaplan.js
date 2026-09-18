import {promiseRequireJS} from '../../utils/promiseRequireJS.js';
import {addTitleLine} from './workers/plotData.js';
import {getScale} from './workers/scale.js';
import {getScatterPositions} from './workers/scatterUtils.js';
import {addDataPointTooltip} from './workers/utils.js';
import {groupAxes} from './components/axis.js';
import {box} from './components/box.js';
import {addClickToEdit, addForm} from './components/form.js';
import {groupLegend} from './components/legend.js';
import {appendClippedPlotAreaG, getSvgDrawPlot, plotInit, plotStart, setTitle} from './components/plotStart.js';
import {scatter} from './components/scatter.js';

const [d3] = await promiseRequireJS(['d3']);

export default function (args) {
    const plotSettings = {
        hasForm: true,
        plotType: 'kaplan',
        analysisType: 'kaplan',
        containerId: args.containerId,
        updateFormId: args.updateFormId,
        plotInnerHeight: parseFloat(args.plotInnerHeight) || 400,
        plotInnerWidth: parseFloat(args.plotInnerWidth) || 400,
        dotSize: args.dotSize || 2,
        fontsizeX: parseFloat(args.fontsizeX) || 14,
        fontsizeY: parseFloat(args.fontsizeY) || 14,
        fontsizeInline: parseFloat(args.fontsizeInline) || 16,
        fontsizeRuler: parseFloat(args.fontsizeRuler) || 12,
        fontsizeTitle: parseFloat(args.fontsizeTitle) || 17,
        fontsizeTitleSub: parseFloat(args.fontsizeTitleSub) || 12,
        titleFormat: args.titleFormat || 'default',
        fontsizeLegend: parseFloat(args.fontsizeLegend) || 12,
        fontsizeLegendHeader: parseFloat(args.fontsizeLegendHeader) || 14,
        axisWidth: parseFloat(args.axisWidth) || 1,
        lineWidth: parseFloat(args.lineWidth) || 1.5,
        addScanPlot: args.addScanPlot === 'yes',
        scanMode: args.scanMode,
        rotateXLabels: args.rotateXLabels ?? 270,
        rotateYLabels: args.rotateYLabels ?? 0,
        plotData: args.plotData,
    };
    plotSettings.plotData.isNumericX = true;
    plotSettings.plotData.isNumericY = true;
    plotSettings.plotData.xDomain = [0, d3.max(plotSettings.plotData.data.map(d => d.xValue))];
    plotSettings.plotData.yDomain = [0, 1];
    plotSettings.plotData.legendFilter = new Set();
    const outerContainerEl = document.getElementById(plotSettings.containerId);
    outerContainerEl.innerText = '';
    plotSettings.containerId = `${plotSettings.containerId}-wrapper-container`;
    const containerEl = document.createElement('div');
    containerEl.id = plotSettings.containerId;
    outerContainerEl.appendChild(containerEl);
    const [formDiv, svg] = plotInit(containerEl);
    plotSettings.plot = d3.select(svg);
    addForm(formDiv, plotSettings);
    plotSettings.drawPlot = getSvgDrawPlot(svg, () => kaplan(plotSettings));
    return plotSettings.drawPlot();
}

async function kaplan(plotSettings) {
    plotStart(plotSettings);
    setTitle(plotSettings);
    const plotFull = plotSettings.plot.append('g');
    plotSettings.xScale = getScale('x', plotSettings);
    plotSettings.yScale = getScale('y', plotSettings);
    groupAxes(plotFull, plotSettings, plotSettings.plotData);
    const plotAreaG = appendClippedPlotAreaG(plotFull, plotSettings.plotInnerWidth, plotSettings.plotInnerHeight, `${plotSettings.containerId}-plot-area-clip`);
    for (const data of plotSettings.plotData.groupData.values()) {
        let nrRemaining = data.length;
        let proportion = 1;
        let e = 0;
        const groupG = plotAreaG.append('g')
            .attr('stroke-width', plotSettings.lineWidth)
            .attr('stroke', data[0].color);

        data.sort((a, b) => a.xValue - b.xValue);

        let lastProportion = 1;
        data.forEach((d, i) => {
            const {xValue, status} = d;
            proportion *= (nrRemaining - status) / nrRemaining;
            nrRemaining--;
            const x1 = data[i - 1] ? data[i - 1].xValue : 0;

            if (status) {
                e++;
                groupG.append('line')
                    .attr('x1', plotSettings.xScale(xValue))
                    .attr('x2', plotSettings.xScale(xValue))
                    .attr('y1', plotSettings.yScale(proportion) + .5)
                    .attr('y2', plotSettings.yScale(lastProportion) - .5)
                    .attr('class', proportion)
                    .datum(d)
                    .call(addDataPointTooltip);
                groupG.append('line')
                    .attr('x1', plotSettings.xScale(x1))
                    .attr('x2', plotSettings.xScale(xValue))
                    .attr('y1', plotSettings.yScale(lastProportion))
                    .attr('y2', plotSettings.yScale(lastProportion));
            } else {
                groupG.append('line')
                    .attr('x1', plotSettings.xScale(xValue))
                    .attr('x2', plotSettings.xScale(xValue))
                    .attr('y1', plotSettings.yScale(proportion) - 3 * plotSettings.lineWidth)
                    .attr('y2', plotSettings.yScale(proportion) + 3 * plotSettings.lineWidth)
                    .attr('class', proportion)
                    .datum(d)
                    .call(addDataPointTooltip);
                groupG.append('line')
                    .attr('x1', plotSettings.xScale(x1))
                    .attr('x2', plotSettings.xScale(xValue))
                    .attr('y1', plotSettings.yScale(proportion))
                    .attr('y2', plotSettings.yScale(proportion));
            }
            lastProportion = proportion;

        });
        const byGroup = d3.group(plotSettings.plotData.data, d => d.colorByValue);
        const legendData = Array.from(byGroup, ([groupId, groupDataPoints]) => ({
            id: groupId,
            label: `${groupId} (n=${groupDataPoints.length})`,
            color: groupDataPoints[0].color,
        }));
        legendData.sort((a, b) => a.label.localeCompare(b.label));
        const colorLegendG = plotSettings.plot.append('g').attr('transform', `translate(${plotSettings.plotInnerWidth + 10},0)`);
        groupLegend(colorLegendG, {
            legendData: legendData,
            legendFilter: plotSettings.plotData.legendFilter,
            plotSettings: plotSettings,
        });
    }

    if (plotSettings.plotData.statistics.hasOwnProperty('inlinePValue')) {

        const inlineText = plotAreaG.append('g');
        inlineText.selectAll("text.inline-text")
            .data(plotSettings.plotData.statistics.inlinePValue)
            .enter()
            .append('text')
            .attr('class', 'inline-text')
            .attr('font-size', `${plotSettings.fontsizeInline}px`)
            .attr('transform', (d, i) => `translate(0, ${plotSettings.fontsizeInline * i})`)
            .text(d => d)
            .call(addClickToEdit, plotSettings.openForm, {
                tabName: 'extra-settings',
                focusSelector: '.inline-text-input',
                titleText: 'click to change font size',
            });

        const dim = inlineText.node().getBBox();
        inlineText.attr('transform', `translate(${plotSettings.plotInnerWidth - dim.width - 5},${plotSettings.plotInnerHeight - dim.height + plotSettings.fontsizeInline})`);
    }
    const sideContainerId = `${plotSettings.containerId}-side-wrapper-container`;
    d3.select(`#${sideContainerId}`).remove();
    if (plotSettings.addScanPlot && plotSettings.scanMode) {
        await addSidePlot(sideContainerId, plotSettings);
    }
}

function addSidePlot(sideContainerId, plotSettings) {
    const plotSettingsScan = Object.assign({}, plotSettings);
    plotSettingsScan.updateFormId = null;
    plotSettingsScan.hasForm = true;
    plotSettingsScan.containerId = sideContainerId;
    const outerContainerEl = document.getElementById(plotSettings.containerId).parentElement;
    const containerEl = document.createElement('div');
    containerEl.id = sideContainerId;
    outerContainerEl.appendChild(containerEl);
    const [formDiv, svg] = plotInit(containerEl);
    plotSettingsScan.plot = d3.select(svg);
    addForm(formDiv, plotSettingsScan);
    plotSettingsScan.drawPlot = getSvgDrawPlot(svg, () => {
        plotSettingsScan.isSidePlot = true;
        if (plotSettings.scanMode === 'scan') {
            addScanPlot(plotSettingsScan);
        } else {
            plotSettingsScan.plotInnerWidth = plotSettings.plotInnerWidth / 2;
            plotSettingsScan.plotInnerHeight = plotSettings.plotInnerHeight / 2;
            addBoxPlot(plotSettingsScan);
        }
    });
    return plotSettingsScan.drawPlot();
}

function addScanPlot(plotSettingsScan) {
    setTitle(plotSettingsScan, {
        titleData: [
            addTitleLine({
                'text': 'Feature values',
                'fontWeight': 600,
                'fontSize': plotSettingsScan.fontsizeTitle,
            }),
        ],
    });
    const expressionDomain = d3.extent(plotSettingsScan.plotData.data.map((d) => d.featureValue));

    plotSettingsScan.plotData.data.sort((a, b) => a.featureValue - b.featureValue);

    const xDomain = plotSettingsScan.plotData.data.map(d => d.id);

    const bestPPoint = plotSettingsScan.plotData.data[plotSettingsScan.plotData.bestPValue.order];

    const width = plotSettingsScan.plotInnerWidth;
    const sideFull = plotSettingsScan.plot.append('g').attr('class', 'plot-area-full');

    const yScale = d3.scaleLinear()
        .range([plotSettingsScan.plotInnerHeight, 0])
        .domain(expressionDomain)
        .nice();
    const yAxis = d3.axisLeft(yScale).tickValues([bestPPoint.featureValue, ...expressionDomain]);

    yAxis.tickSize(plotSettingsScan.axisWidth * 2).tickFormat(d3.format(",.2r"));
    const yAxisG = sideFull.append("g")
        .attr('class', 'side-y-axis')
        .style('font-size', plotSettingsScan.fontsizeRuler)
        .call(yAxis);

    const bBox = yAxisG.node().getBBox();
    yAxisG.append('text')
        .text('value')
        .attr('fill', '#111111')
        .attr('text-anchor', 'middle')
        .attr('font-size', plotSettingsScan.fontsizeY)
        .attr('font-weight', 600)
        .attr('transform', `translate(${-plotSettingsScan.fontsizeY - bBox.width}, ${plotSettingsScan.plotInnerHeight / 2}) rotate(-90)`);

    sideFull.selectAll('.domain').style('stroke-width', `${plotSettingsScan.axisWidth}px`);

    const xScale = d3.scalePoint()
        .range([0, width])
        .domain(xDomain)
        .padding(0);

    const xAxis = d3.axisBottom(xScale).tickSize(0).tickFormat(() => '');
    sideFull.append("g")
        .attr('class', 'side-x-axis')
        .call(xAxis)
        .attr('transform', `translate(0,${plotSettingsScan.plotInnerHeight})`);

    const plotArea = sideFull.append('g').attr('class', 'plot-area');

    const pData = [];
    plotSettingsScan.plotData.pValues.forEach((d) => {
        const data = plotSettingsScan.plotData.data[d.order];
        if (!data?.tooltipDataUpdated) {
            data.tooltipData['hover_text'] = `p-value: ${d.pValue}</br>value: ${data.featureValue}</br>${data.tooltipData['hover_text']}`;
            data.tooltipDataUpdated = true;
        }
        pData.push({
            'logP': -1 * Math.log(d.pValue) / Math.log(10),
            'pValue': d.pValue,
            'order': d.order,
            ...data,
        });
    });

    plotArea.append('text')
        .text('sorted by feature value')
        .attr('text-anchor', 'middle')
        .attr('transform', `translate(${width / 2},${plotSettingsScan.plotInnerHeight})`)
        .style('dominant-baseline', 'hanging');

    plotArea.append('line')
        .attr('x1', xScale(plotSettingsScan.plotData.data[0].id))
        .attr('x2', xScale(plotSettingsScan.plotData.data.slice(-1)[0].id))
        .attr('y1', yScale(bestPPoint.featureValue))
        .attr('y2', yScale(bestPPoint.featureValue))
        .attr('stroke', '#b3b3b3');

    plotArea.append('line')
        .attr('x1', xScale(bestPPoint.id))
        .attr('x2', xScale(bestPPoint.id))
        .attr('y1', 0)
        .attr('y2', plotSettingsScan.plotInnerHeight)
        .attr('stroke', '#b3b3b3');

    const updatePValueCutoff = (order) => {
        const el = document.querySelector('[name=selected_a]');
        el.value = order;
        const form = document.querySelector('form[name=rekaplan]');
        form.submit();
    };

    plotArea.selectAll('g.point').data(plotSettingsScan.plotData.data).enter()
        .append('g')
        .attr('class', 'point')
        .append('circle')
        .attr('r', plotSettingsScan.dotSize)
        .attr('cx', (d) => xScale(d.id))
        .attr('cy', (d) => yScale(d.featureValue))
        .attr('fill', (d) => d.status ? '#c80000' : '#00c800')
        .on('click', (e, d) => {
            const dp = pData.find((pd) => pd.id === d.id);
            updatePValueCutoff(dp.order);
        })
        .call(addDataPointTooltip);

    const pPlotG = sideFull.append('g').attr('transform', `translate(0,${plotSettingsScan.plotInnerHeight + 10})`);
    const domainMax = d3.max(pData.map(d => d.logP));
    const pYScale = d3.scaleLinear().range([50, 0]).domain([-0, domainMax]);
    const pYAxis = d3.axisLeft(pYScale).tickValues([0, domainMax]);
    const pYAxisG = pPlotG.append('g')
        .attr('class', 'side-p-y-axis')
        .style('font-size', plotSettingsScan.fontsizeRuler)
        .call(pYAxis);

    pYAxisG.append('text')
        .text('-log p-value')
        .attr('fill', '#111111')
        .attr('text-anchor', 'middle')
        .attr('font-size', "10px")
        .attr('font-weight', 400)
        .attr('transform', `translate(${-plotSettingsScan.fontsizeY - bBox.width}, 25) rotate(-90)`);

    pPlotG.append('g')
        .attr('class', 'side-p-x-axis')
        .style('font-size', plotSettingsScan.fontsizeRuler).call(xAxis)
        .attr('transform', 'translate(0,50)');

    const plotAreaP = pPlotG.append('g');

    plotAreaP.selectAll('g.p-rect').data(pData).enter()
        .append('g')
        .attr('class', 'p-rect')
        .append('line')
        .attr('x1', (d) => xScale(d.id))
        .attr('x2', (d) => xScale(d.id))
        .attr('y2', (d) => pYScale(d.logP))
        .attr('y1', pYScale(-0))
        .attr('stroke', '#7393B3')
        .attr('stroke-width', plotSettingsScan.dotSize * 2)
        .on('click', (e, d) => updatePValueCutoff(d.order))
        .call(addDataPointTooltip);
}


function addBoxPlot(plotSettingsScan) {
    plotSettingsScan.defaultSurfaceColor = '#cccccc';
    plotSettingsScan.boxStrokeWidth = 1;
    plotSettingsScan.boxOpacity = .7;
    plotSettingsScan.datapointOpacity = .7;
    plotSettingsScan.boxWidth = 75;
    plotSettingsScan.groupsOrderBy = 'average';
    plotSettingsScan.plotData.isNumericX = false;
    plotSettingsScan.altYValues = 'featureValue';
    plotSettingsScan.colorMode = 'givenColor';
    plotSettingsScan.plotData.xLabel = 'groups';
    plotSettingsScan.plotData.yLabel = 'feature values';
    plotSettingsScan.plotData.xDomain = plotSettingsScan.plotData.sumStats.map(d => d.key);
    plotSettingsScan.plotData.yDomain = d3.extent(plotSettingsScan.plotData.data.map(d => d.featureValue));

    plotStart(plotSettingsScan);
    setTitle(plotSettingsScan, {
        titleData: [
            addTitleLine({
                'text': 'Feature values per group',
                'fontWeight': 600,
                'fontSize': plotSettingsScan.fontsizeTitle,
            }),
        ],
    });
    const plotFull = plotSettingsScan.plot.append('g').attr('class', 'plot-area-full');
    plotSettingsScan.xScale = getScale('x', plotSettingsScan);
    plotSettingsScan.yScale = getScale('y', plotSettingsScan);
    groupAxes(plotFull, plotSettingsScan, plotSettingsScan.plotData);
    appendClippedPlotAreaG(plotFull, plotSettingsScan.plotInnerWidth, plotSettingsScan.plotInnerHeight, `${plotSettingsScan.containerId}-plot-area-clip`);
    plotFull.call(box, plotSettingsScan);
    getScatterPositions(plotSettingsScan);
    scatter(plotSettingsScan);
}
