import { useColorMode } from '@docusaurus/theme-common';
import React, { useEffect, useState } from 'react';
import useBaseUrl from '@docusaurus/useBaseUrl';
import Plot from 'react-plotly.js';
import Parser from 'papaparse';

function CostPlot() {
  const { colorMode, setColorMode } = useColorMode();

  const [data, setData] = useState([]);
  const filename = useBaseUrl('/data/cost.csv')
  const singleSampleColor = colorMode == 'dark' ? '#25c2a0' : '#0f4e40'
  const multiSampleColor = colorMode == 'dark' ? '#794fef' : '#260b75'

  // Parse the data asynchronously.
  useEffect(() => {
    fetch(filename)
      .then((response) => response.text())
      .then((data) => {
        const parsedData = Parser.parse(data, {
          header: false,
          skipEmptyLines: true,
          fastMode: true,
          delimiter: ','
        }).data;

        const rows = parsedData;

        const plotData = [];
        //plotData.push({y: [0], name: 'Single-sample pipeline', marker: {color: singleSampleColor}})
        //plotData.push({y: [0], name: 'Multi-sample pipeline', marker: {color: multiSampleColor}})

        for (let i = 0; i < rows.length; i++) {
          const [groupName, seriesName, ...rowValues] = rows[i];
          const y = rowValues.map(Number); // Convert string values to numbers

          plotData.push({
            x: seriesName, // Array(y.length).fill(groupName),
            y: y,
            type: 'box',
            boxpoints: 'all',
            jitter: 0.8,
            whiskerwidth: 0.5,
            //showlegend: false,
            name: groupName + " - " + seriesName,
            marker: {
              size: 5,
              color: groupName == 'Single-sample' ? singleSampleColor : multiSampleColor,
            },
            line: {
              width: 2,
            }
          });
        }

        setData(plotData);
      });
  }, []);

  if (data)
    return (
     <div style={{ display: 'flex', justifyContent: 'center' }}>
      <Plot
        layout={{
          title: 'Workflow Execution Cost',
          font: {color: colorMode == 'dark' ? 'rgba(255,255,255,0.5)' : 'rgb(0,0,0)'},
          xaxis: {title: 'Cohort'},
          yaxis: {
            title: 'Cost (USD)',
            gridcolor: colorMode == 'dark' ? 'rgba(200,200,200,0.1)' : 'rgba(255,255,255,1)'
          },

          autosize: true,
          plot_bgcolor: colorMode == 'dark' ? 'rgba(200,200,200,0.1)' : 'rgb(234,234,234)',
          paper_bgcolor: colorMode == 'dark' ? 'rgba(234,234,234,0.0)' : 'rgb(255,255,255)',
          //boxmode: 'group',
          //legend: {orientation: 'h', xanchor: 'center', yanchor: 'center', x:0.5, y:-0.5},
          margin: {t:50, l:50, r:50, b:50}
        }}
        useResizeHandler={true}
        style={{ width: '100%', height: '100%' }}
        data={data}
      />
      </div>
    );
}

export default CostPlot;