import {fetchGetJson, fetchPost} from './utils/fetch.js';
import {showDatasetCard} from './general/datasetCard.js';
import {promiseRequireJS} from './utils/promiseRequireJS.js';

const [$] = await promiseRequireJS(['jquery', 'jqx-all']);

export default function ({inputId, filterDatasets, multiple}) {
    const datasetSelect = document.getElementById(inputId);
    $(datasetSelect).on('mousedown', event => event.preventDefault());
    $(datasetSelect).on('click', event => {
        event.preventDefault();
        let $loaderDiv = $("<div>")
            .appendTo(document.body)
            .jqxLoader({isModal: true, imagePosition: 'center', autoOpen: true});
        fetchGetJson({
            'json_option': 'json_dataset',
            'selected': [...datasetSelect.selectedOptions].map(option => option.value).join(';'),
            'filter_datasets': filterDatasets ? filterDatasets.join(';') : '',
        }).then(d => {
            $loaderDiv.jqxLoader('close');
            openDatasetChanger(d);
        });
    });

    function openDatasetChanger(dataRows) {
        const $modalDiv = $("<div>").appendTo(document.body);
        let currentSelectionText = 'None';
        if (datasetSelect.multiple) {
            currentSelectionText = datasetSelect.selectedOptions.length;
        } else if (datasetSelect.selectedOptions[0].value !== '') {
            currentSelectionText = datasetSelect.selectedOptions[0].textContent;
        }
        $modalDiv.append(`<div>Data set selection - current: <b>${currentSelectionText}</b></div>`);
        const $modalContentDiv = $("<div>").appendTo($modalDiv);
        $modalDiv.jqxWindow({
            minWidth: "98%",
            height: "98%",
            resizable: false,
            isModal: true,
        });
        $modalDiv.on('close', () => $modalDiv.jqxWindow('destroy'));
        const $gridDiv = $("<div>").appendTo($modalContentDiv).css({'min-height': '400px', 'margin': '20px'});
        $gridDiv.jqxGrid({
            width: '97%',
            height: '50%',
            columnsResize: true,
            sortable: true,
            showFilterRow: true,
            filterRowHeight: 32,
            filterable: true,
            autoShowFilterIcon: false,
            selectionMode: multiple ? 'checkbox' : 'singlerow',
            editable: true,
            showstatusbar: true,
            renderstatusbar: statusbar => {
                const $bottomBarDiv = $("<div>").appendTo(statusbar).css({
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'space-between',
                    'height': '100%',
                });
                $("<input type='button' value='Confirm selection'>").appendTo($bottomBarDiv).css({
                    'margin-left': '10px',
                    'height': '30px',
                    'width': '150px',
                    'font-weight': 'bold',
                    'color': '#000088',
                    'background-color': '#c8daf6',
                    'border': '2px solid #2779FF',
                    'border-radius': '5px',
                }).on('click', () => {
                    datasetSelect.textContent = '';
                    const selectedRowIndexes = $gridDiv.jqxGrid('getselectedrowindexes');
                    if (selectedRowIndexes.length) {
                        datasetSelect.multiple = !!multiple;
                        selectedRowIndexes.forEach(selectedRowIndex => {
                            const rowData = $gridDiv.jqxGrid('getrowdata', selectedRowIndex);
                            datasetSelect.add(new Option(rowData['showname'], rowData['dataset'], true, true));
                        });
                        datasetSelect.size = Math.min(selectedRowIndexes.length, 10);
                    } else {
                        datasetSelect.multiple = false;
                        const optionText = multiple ? 'Select data sets' : 'Select a data set';
                        datasetSelect.add(new Option(optionText, '', true, true));
                        datasetSelect.size = 1;
                    }
                    datasetSelect.dispatchEvent(new Event('change'));
                    $modalDiv.jqxWindow('destroy');
                });
                const $rowCountSpan = $("<span>").appendTo($bottomBarDiv).text(`rows: ${dataRows.length}`).css({'margin-right': '10px'});
                $gridDiv.on('filter', () => $rowCountSpan.text(`rows: ${$gridDiv.jqxGrid('getrows').length}`));
            },
            source: new $.jqx.dataAdapter({
                dataFields: [
                    {name: 'dataset_samples', type: 'float'},
                    {name: 'species', type: 'string'},
                    {name: 'dataset', type: 'string'},
                    {name: 'showname', type: 'string'},
                    {name: 'datatype', type: 'string'},
                    {name: 'dataset_id', type: 'string'},
                    {name: 'dataset_class', type: 'string'},
                    {name: 'dataset_composition', type: 'string'},
                    {name: 'dataset_material', type: 'string'},
                    {name: 'dataset_subclasses', type: 'string'},
                    {name: 'dataset_author', type: 'string'},
                    {name: 'dataset_normalization', type: 'string'},
                    {name: 'platform', type: 'string'},
                    {name: 'dataset_r2_date', type: 'date'},
                    {name: 'dataset_date', type: 'date'},
                    {name: 'access', type: 'string'},
                    {name: 'favourite', type: 'bool'},
                ],
                id: 'dataset',
                localData: dataRows,
            }),
            columns: [
                {
                    text: 'Species',
                    dataField: 'species',
                    width: 70,
                    filterType: 'checkedlist',
                    editable: false,
                },
                {
                    text: 'Data type',
                    dataField: 'datatype',
                    width: 100,
                    filterType: 'checkedlist',
                    editable: false,
                },
                {
                    text: 'Category',
                    dataField: 'dataset_class',
                    width: 100,
                    filterType: 'checkedlist',
                    editable: false,
                },
                {
                    text: 'Tissue/Tumor',
                    dataField: 'dataset_subclasses',
                    editable: false,
                },
                {
                    text: 'Author',
                    dataField: 'dataset_author',
                    width: 120,
                    editable: false,
                },
                {
                    text: 'N',
                    dataField: 'dataset_samples',
                    filterType: 'number',
                    width: 70,
                    editable: false,
                },
                {
                    text: 'Normalization',
                    dataField: 'dataset_normalization',
                    filterType: 'checkedlist',
                    width: 100,
                    editable: false,
                },
                {
                    text: 'Platform',
                    dataField: 'platform',
                    width: 100,
                    editable: false,
                },
                {
                    text: 'Composition',
                    dataField: 'dataset_composition',
                    width: 65,
                    filterType: 'checkedlist',
                    editable: false,
                },
                {
                    text: 'Material',
                    dataField: 'dataset_material',
                    width: 65,
                    filterType: 'checkedlist',
                    editable: false,
                },
                {
                    text: 'Accession',
                    dataField: 'dataset_id',
                    width: 75,
                    editable: false,
                },
                {
                    text: 'Release date',
                    dataField: 'dataset_date',
                    filterType: 'date',
                    width: 90,
                    cellsFormat: 'yyyy-MM-dd',
                    editable: false,
                },
                {
                    text: 'R2 date',
                    dataField: 'dataset_r2_date',
                    filterType: 'date',
                    width: 90,
                    cellsFormat: 'yyyy-MM-dd',
                    editable: false,
                },
                {
                    text: 'Access',
                    dataField: 'access',
                    width: 70,
                    filterType: 'checkedlist',
                    editable: false,
                },
                {
                    text: 'Favourite',
                    dataField: 'favourite',
                    columnType: 'checkbox',
                    filterType: 'bool',
                    width: 57,
                    editable: true,
                    cellClassName: (row, column, value) => value ? 'greenCell' : '',
                },
            ],
        });
        $modalContentDiv.append("<div id='info'></div>");
        if (datasetSelect.selectedOptions[0].value !== '') {
            const selectedRowIndexes = [...datasetSelect.selectedOptions].map(option => $gridDiv.jqxGrid('getrowboundindexbyid', option.value)).sort();
            selectedRowIndexes.forEach(rowIndex => $gridDiv.jqxGrid('selectrow', rowIndex));
            showDatasetCard({containerId: 'info', dataset: $gridDiv.jqxGrid('getrowid', selectedRowIndexes[0])});
        }
        $gridDiv.on('rowselect', event => {
            const rowData = event.args.row;
            if (rowData) {
                showDatasetCard({containerId: 'info', dataset: rowData['dataset']});
            }
        });
        if (!('favourite' in dataRows[0])) {
            $gridDiv.jqxGrid('hideColumn', 'favourite');
        } else {
            $gridDiv.on('cellvaluechanged', event => {
                const args = event.args;
                if (args.datafield === 'favourite') {
                    fetchPost({
                        json_option: 'json_toggle_favourite_dataset',
                        dataset: $gridDiv.jqxGrid('getRowData', args.rowindex)['dataset'],
                        toggle: args.newvalue,
                    });
                }
            });
        }
    }
}

