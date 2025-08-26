<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Enhanced Protein-Ligand Contact Analysis</title>
    <style>
        body {
            margin: 0; 
            padding: 20px; 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }
        .container {
            max-width: 1400px; 
            margin: 0 auto; 
            background: rgba(255, 255, 255, 0.98); 
            border-radius: 20px; 
            box-shadow: 0 20px 60px rgba(0,0,0,0.3), 0 0 100px rgba(102, 126, 234, 0.1);
            padding: 35px;
            backdrop-filter: blur(10px);
        }
        h1 {
            text-align: center; 
            color: #2C3E50; 
            font-size: 32px; 
            margin-bottom: 10px; 
            font-weight: 700;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            text-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .subtitle {
            text-align: center; 
            color: #7F8C8D; 
            font-size: 16px; 
            margin-bottom: 30px;
            font-weight: 300;
            letter-spacing: 0.5px;
        }
        #canvas {
            border: 1px solid rgba(102, 126, 234, 0.2);
            border-radius: 15px; 
            background: linear-gradient(to bottom, #ffffff, #f8f9fa);
            cursor: grab; 
            display: block; 
            margin: 0 auto;
            box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        }
        #canvas:active { cursor: grabbing; }
        .controls {
            display: flex; 
            justify-content: center; 
            gap: 12px; 
            margin-bottom: 20px; 
            flex-wrap: wrap;
            align-items: center;
        }
        .control-btn {
            padding: 12px 24px; 
            border: none; 
            border-radius: 10px; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; 
            font-size: 14px; 
            cursor: pointer; 
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            font-weight: 500;
            letter-spacing: 0.3px;
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);
        }
        .control-btn:hover {
            transform: translateY(-2px) scale(1.02);
            box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4);
        }
        .control-btn.secondary {
            background: linear-gradient(135deg, #95a5a6 0%, #7f8c8d 100%);
            box-shadow: 0 4px 15px rgba(149, 165, 166, 0.3);
        }
        .control-btn.secondary:hover {
            box-shadow: 0 6px 20px rgba(149, 165, 166, 0.4);
        }
        .control-btn.active {
            background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%);
            box-shadow: 0 4px 15px rgba(231, 76, 60, 0.3);
        }
        .zoom-controls {
            display: flex; 
            justify-content: center; 
            gap: 10px; 
            margin-bottom: 15px; 
        }
        .zoom-btn {
            width: 45px; 
            height: 45px; 
            border-radius: 50%; 
            border: 2px solid transparent;
            background: linear-gradient(white, white) padding-box,
                        linear-gradient(135deg, #667eea 0%, #764ba2 100%) border-box;
            color: #667eea; 
            font-size: 20px; 
            cursor: pointer; 
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            font-weight: 600;
            box-shadow: 0 2px 10px rgba(102, 126, 234, 0.2);
        }
        .zoom-btn:hover {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) padding-box,
                        linear-gradient(135deg, #667eea 0%, #764ba2 100%) border-box;
            color: white;
            transform: scale(1.1);
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
        }
        .info-panel {
            background: linear-gradient(135deg, #f5f7fa 0%, #ecf0f1 100%);
            border-radius: 15px; 
            padding: 25px; 
            margin-top: 25px;
            box-shadow: inset 0 2px 10px rgba(0,0,0,0.05);
        }
        .file-info {
            background: linear-gradient(135deg, #dfe6e9 0%, #b2bec3 100%);
            border-radius: 10px; 
            padding: 15px; 
            margin-bottom: 15px; 
            font-size: 13px;
            color: #2c3e50;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }
        .stats {
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); 
            gap: 15px; 
            margin: 20px 0; 
        }
        .stat-item {
            background: white; 
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
            transition: all 0.3s ease;
            border: 1px solid rgba(102, 126, 234, 0.1);
        }
        .stat-item:hover {
            transform: translateY(-3px);
            box-shadow: 0 8px 25px rgba(102, 126, 234, 0.15);
        }
        .stat-value {
            font-size: 28px; 
            font-weight: bold; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .stat-label {
            font-size: 12px; 
            color: #7F8C8D; 
            margin-top: 8px;
            text-transform: uppercase;
            letter-spacing: 1px;
        }
        .legend {
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); 
            gap: 15px; 
            margin-top: 20px; 
        }
        .legend-item {
            display: flex; 
            align-items: center; 
            gap: 12px; 
            font-size: 14px; 
            color: #2C3E50;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s;
        }
        .legend-item:hover {
            background: rgba(102, 126, 234, 0.05);
        }
        .legend-color {
            width: 24px; 
            height: 24px; 
            border-radius: 50%;
            box-shadow: 0 2px 5px rgba(0,0,0,0.2);
        }
        .color-bar {
            margin: 20px 0;
            padding: 20px;
            background: white;
            border-radius: 12px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.08);
        }
        .color-gradient-container {
            position: relative;
            width: 100%;
            margin: 15px 0 25px 0;  /* 底部增加空间 */
        }
        .color-gradient {
            height: 40px;
            width: 100%;
            background: linear-gradient(to right, 
                #95a5a6 0%, 
                #3498db 25%, 
                #9b59b6 50%, 
                #e74c3c 75%, 
                #c0392b 100%);
            border-radius: 20px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.15), inset 0 2px 4px rgba(255,255,255,0.3);
        }
        .color-labels {
            display: flex;
            justify-content: space-between;
            width: 100%;
            margin-top: 8px;
            font-size: 12px;
            color: #7F8C8D;
            font-weight: 600;
        }
        .color-bar-title {
            font-size: 16px;
            font-weight: bold;
            color: #2C3E50;
            margin-bottom: 10px;
        }
        .tooltip {
            position: absolute; 
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
            color: white; 
            padding: 12px 16px; 
            border-radius: 8px; 
            font-size: 13px; 
            pointer-events: none; 
            z-index: 1000; 
            display: none;
            box-shadow: 0 8px 25px rgba(0,0,0,0.3);
            border: 1px solid rgba(255,255,255,0.1);
        }
        .help-text {
            text-align: center; 
            color: #7F8C8D; 
            font-size: 13px; 
            margin: 15px 0;
            padding: 10px;
            background: rgba(102, 126, 234, 0.05);
            border-radius: 8px;
        }
        .export-controls {
            display: inline-flex;
            gap: 10px;
            align-items: center;
            background: rgba(255,255,255,0.5);
            padding: 8px 12px;
            border-radius: 10px;
        }
        .export-controls select, .export-controls input {
            color: #000000; 
            padding: 10px;
            border: 1px solid rgba(102, 126, 234, 0.3);
            border-radius: 8px;
            font-size: 14px;
            transition: all 0.3s;
        }
        .export-controls select:focus, .export-controls input:focus {
            outline: none;
            border-color: #667eea;
            box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1);
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>Enhanced Protein-Ligand Contact Analysis</h1>
        <p class="subtitle">Interactive visualization with improved 3D mode and alignment</p>
        
        <div class="file-info">
            <strong>Analysis Details:</strong> Trajectory: md.xtc | Topology: md.gro | Ligand: 3q4j_ligand.sdf (LIG) | Frames: 1000
        </div>
        
        <div class="zoom-controls">
            <button class="zoom-btn" onclick="zoomIn()">+</button>
            <button class="zoom-btn" onclick="zoomOut()">−</button>
            <button class="zoom-btn" onclick="centerView()" style="width: auto; padding: 0 15px;">Center</button>
        </div>
        
        <div class="controls">
            <button class="control-btn" onclick="resetPositions()">Reset Layout</button>
            <button class="control-btn" onclick="rotateWheel()">Rotate Wheel</button>
            <button class="control-btn" onclick="toggleDistanceLock()" id="distanceLockBtn">🔒 Distance Locked</button>
            <button class="control-btn" onclick="toggle3DMode()" id="mode3DBtn">2D Mode</button>
            <button class="control-btn secondary" onclick="toggleConnections()" id="connectBtn">Hide Connections</button>
            <button class="control-btn secondary" onclick="toggleHydrogens()" id="hydrogenBtn">Hide H atoms</button> 
            <div class="export-controls" style="display: inline-flex; gap: 10px; align-items: center;">
                <select id="exportQuality" class="control-btn secondary" style="padding: 10px;">
                    <option value="1">Standard (1x)</option>
                    <option value="2">High (2x)</option>
                    <option value="4">Ultra HD (4x)</option>
                    <option value="8">8K (8x)</option>
                </select>
                <input type="text" id="exportFilename" placeholder="filename" value="contact_analysis" style="padding: 8px; border: 1px solid #95A5A6; border-radius: 5px;">
                <button class="control-btn secondary" onclick="exportImage()">Export PNG</button>
            </div>
        </div>

        <p class="help-text">🖱️ In 2D: Drag residues to rotate | In 3D: Drag to rotate molecule | Middle-click/Shift+Drag to pan | Scroll to zoom</p>

        <canvas id="canvas" width="1200" height="800"></canvas>
        
        <div class="info-panel">
            <h3 style="margin-top: 0;">Contact Analysis Results</h3>
            <div class="stats">
                <div class="stat-item"><div class="stat-value">15</div><div class="stat-label">Total Contacts</div></div>
                <div class="stat-item"><div class="stat-value">2</div><div class="stat-label">High Frequency</div></div>
                <div class="stat-item"><div class="stat-value">71%</div><div class="stat-label">Max Frequency</div></div>
                <div class="stat-item"><div class="stat-value">96</div><div class="stat-label">Ligand Atoms</div></div>
            </div>
            
            <div class="color-bar">
                <div class="color-bar-title">Contact Frequency</div>
                <div class="color-gradient-container">
                    <div class="color-gradient"></div>
                </div>
                <div class="color-labels">
                    <span>0%</span>
                    <span>25%</span>
                    <span>50%</span>
                    <span>75%</span>
                    <span>100%</span>
                </div>
            </div>
            <div class="legend">
                <div class="legend-item"><div class="legend-color" style="background: #27ae60;"></div><span>Carbon atoms</span></div>
                <div class="legend-item"><div class="legend-color" style="background: #e91e63;"></div><span>Oxygen atoms</span></div>
                <div class="legend-item"><div class="legend-color" style="background: #2196f3;"></div><span>Nitrogen atoms</span></div>
                <div class="legend-item"><div class="legend-color" style="background: #ffc107;"></div><span>Sulfur atoms</span></div>
            </div>
        </div>
    </div>

    <div class="tooltip" id="tooltip"></div>

    <script>
        const LIGAND_ATOMS = [
  {
    "id": "L0",
    "x": -264.98041721013004,
    "y": 13.724346761670585,
    "x3d": -40.792172062884774,
    "y3d": 9.986192058553929,
    "z3d": -9.062913828029673,
    "element": "C",
    "radius": 18,
    "contacts": 174,
    "trajIndex": 11372
  },
  {
    "id": "L1",
    "x": -281.0944441901936,
    "y": -14.850225142135551,
    "x3d": -40.16403665693611,
    "y3d": 14.75225300097798,
    "z3d": -11.27345019273257,
    "element": "O",
    "radius": 18,
    "contacts": 2522,
    "trajIndex": 11373
  },
  {
    "id": "L2",
    "x": -292.42552129373684,
    "y": 49.386178765088346,
    "x3d": -46.60093210291175,
    "y3d": 7.14142469168948,
    "z3d": -8.437766617962856,
    "element": "C",
    "radius": 18,
    "contacts": 70,
    "trajIndex": 11374
  },
  {
    "id": "L3",
    "x": -220.373812707874,
    "y": 19.66158810587317,
    "x3d": -36.33472985390046,
    "y3d": 6.730653674853867,
    "z3d": -6.92775191388982,
    "element": "N",
    "radius": 18,
    "contacts": 106,
    "trajIndex": 11375
  },
  {
    "id": "L4",
    "x": -192.92870862426724,
    "y": -16.000243897544596,
    "x3d": -30.28913417617187,
    "y3d": 8.480656492281813,
    "z3d": -7.405989953521354,
    "element": "C",
    "radius": 18,
    "contacts": 31,
    "trajIndex": 11376
  },
  {
    "id": "L5",
    "x": -157.2668766208495,
    "y": 11.444860186062158,
    "x3d": -26.518862777739066,
    "y3d": 3.0827934211550714,
    "z3d": -6.438103318411846,
    "element": "C",
    "radius": 18,
    "contacts": 18,
    "trajIndex": 11377
  },
  {
    "id": "L6",
    "x": -165.52669261145402,
    "y": 43.19298556066263,
    "x3d": -28.563876094394203,
    "y3d": -1.4624084181730992,
    "z3d": -4.582515472786026,
    "element": "O",
    "radius": 18,
    "contacts": 239,
    "trajIndex": 11378
  },
  {
    "id": "L7",
    "x": -228.59054062768502,
    "y": -43.44534798115135,
    "x3d": -28.733539796734075,
    "y3d": 13.51070082969297,
    "z3d": -3.3420850709563834,
    "element": "C",
    "radius": 18,
    "contacts": 83,
    "trajIndex": 11379
  },
  {
    "id": "L8",
    "x": -264.2523726311028,
    "y": -70.89045206475812,
    "x3d": -30.317142280358638,
    "y3d": 12.53913386564301,
    "z3d": 3.009831162936019,
    "element": "C",
    "radius": 18,
    "contacts": 726,
    "trajIndex": 11380
  },
  {
    "id": "L9",
    "x": -299.91420463452056,
    "y": -98.33555614836489,
    "x3d": -26.117697958817942,
    "y3d": 8.733764841222433,
    "z3d": 6.274901065297654,
    "element": "C",
    "radius": 18,
    "contacts": 158,
    "trajIndex": 11381
  },
  {
    "id": "L10",
    "x": -341.51327798214095,
    "y": -81.1740557297155,
    "x3d": -21.151140872203232,
    "y3d": 7.729527282848438,
    "z3d": 4.7604385499286614,
    "element": "O",
    "radius": 18,
    "contacts": 672,
    "trajIndex": 11382
  },
  {
    "id": "L11",
    "x": -293.976963290318,
    "y": -142.94216065062093,
    "x3d": -28.237906060549168,
    "y3d": 6.607660196989375,
    "z3d": 11.343121701720555,
    "element": "N",
    "radius": 18,
    "contacts": 1952,
    "trajIndex": 11383
  },
  {
    "id": "L12",
    "x": -115.66780327322913,
    "y": -5.716640232587142,
    "x3d": -20.79812205299092,
    "y3d": 3.765812203056614,
    "z3d": -7.827536291365697,
    "element": "N",
    "radius": 18,
    "contacts": 915,
    "trajIndex": 11384
  },
  {
    "id": "L13",
    "x": -80.00597126981135,
    "y": 21.72846385101962,
    "x3d": -16.424521146133284,
    "y3d": -0.4269997447474065,
    "z3d": -6.0373961843618495,
    "element": "C",
    "radius": 18,
    "contacts": 24,
    "trajIndex": 11385
  },
  {
    "id": "L14",
    "x": -52.5608671862046,
    "y": -13.933368152398145,
    "x3d": -10.937584400415352,
    "y3d": 3.137470475892428,
    "z3d": -4.7475040763336915,
    "element": "C",
    "radius": 18,
    "contacts": 234,
    "trajIndex": 11386
  },
  {
    "id": "L15",
    "x": -69.72236760485391,
    "y": -55.532441500018514,
    "x3d": -10.072555869417727,
    "y3d": 7.878895214738926,
    "z3d": -6.981874548576787,
    "element": "O",
    "radius": 18,
    "contacts": 1024,
    "trajIndex": 11387
  },
  {
    "id": "L16",
    "x": -107.4510753534181,
    "y": 57.39029585443739,
    "x3d": -15.46014777402609,
    "y3d": -4.893069846403727,
    "z3d": -10.855680256344883,
    "element": "C",
    "radius": 18,
    "contacts": 587,
    "trajIndex": 11388
  },
  {
    "id": "L17",
    "x": -134.89617943702487,
    "y": 93.05212785785517,
    "x3d": -11.412803871973079,
    "y3d": -9.98059146458062,
    "z3d": -9.405640975308033,
    "element": "C",
    "radius": 18,
    "contacts": 153,
    "trajIndex": 11389
  },
  {
    "id": "L18",
    "x": -162.3412835206316,
    "y": 128.7139598612729,
    "x3d": -10.618692023232096,
    "y3d": -13.817649745533252,
    "z3d": -14.747653684887668,
    "element": "C",
    "radius": 18,
    "contacts": 350,
    "trajIndex": 11390
  },
  {
    "id": "L19",
    "x": -99.23434743360711,
    "y": 120.4972319414619,
    "x3d": -13.648617297572772,
    "y3d": -13.648855246129363,
    "z3d": -4.38484361354521,
    "element": "C",
    "radius": 18,
    "contacts": 177,
    "trajIndex": 11391
  },
  {
    "id": "L20",
    "x": -7.9542626839485315,
    "y": -7.996126808195536,
    "x3d": -7.064953760179568,
    "y3d": 0.5951675808499286,
    "z3d": -1.059536540199889,
    "element": "N",
    "radius": 18,
    "contacts": 66,
    "trajIndex": 11392
  },
  {
    "id": "L21",
    "x": 19.49084139965823,
    "y": -43.657958811613305,
    "x3d": -1.5127976014956883,
    "y3d": 3.2264894450985624,
    "z3d": 0.3898207476674662,
    "element": "C",
    "radius": 18,
    "contacts": 560,
    "trajIndex": 11393
  },
  {
    "id": "L22",
    "x": 55.152673403075994,
    "y": -16.212854728006544,
    "x3d": 2.7226274196611695,
    "y3d": -1.8365585720327284,
    "z3d": 1.2600531806533062,
    "element": "C",
    "radius": 18,
    "contacts": 142,
    "trajIndex": 11394
  },
  {
    "id": "L23",
    "x": 43.822335967404555,
    "y": 27.337385840718472,
    "x3d": 1.025992003108139,
    "y3d": -6.784180939003835,
    "z3d": 2.2575663823738665,
    "element": "O",
    "radius": 18,
    "contacts": 312,
    "trajIndex": 11395
  },
  {
    "id": "L24",
    "x": -16.170990603759535,
    "y": -71.10306289522006,
    "x3d": -2.0110151414982655,
    "y3d": 6.826204104312642,
    "z3d": 5.891256729159864,
    "element": "C",
    "radius": 18,
    "contacts": 37,
    "trajIndex": 11396
  },
  {
    "id": "L25",
    "x": -51.832822607177306,
    "y": -98.54816697882684,
    "x3d": 2.9820013044032088,
    "y3d": 10.918365373487111,
    "z3d": 6.384738879566501,
    "element": "C",
    "radius": 18,
    "contacts": 591,
    "trajIndex": 11397
  },
  {
    "id": "L26",
    "x": -93.43189595479768,
    "y": -81.3866665601775,
    "x3d": 6.875180703619237,
    "y3d": 11.57671502709111,
    "z3d": 2.8978595698504352,
    "element": "O",
    "radius": 18,
    "contacts": 1400,
    "trajIndex": 11398
  },
  {
    "id": "L27",
    "x": -45.895581262974694,
    "y": -143.1547714810829,
    "x3d": 2.878024437108591,
    "y3d": 13.957742868383729,
    "z3d": 11.344925493081607,
    "element": "O",
    "radius": 18,
    "contacts": 1303,
    "trajIndex": 11399
  },
  {
    "id": "L28",
    "x": 96.75174675069637,
    "y": -33.37435514665585,
    "x3d": 8.466756889993949,
    "y3d": -0.4412405617083337,
    "z3d": 0.9192209527684091,
    "element": "N",
    "radius": 18,
    "contacts": 1009,
    "trajIndex": 11400
  },
  {
    "id": "L29",
    "x": 132.41357875411413,
    "y": -5.92925106304909,
    "x3d": 13.013306476094666,
    "y3d": -4.828181145374453,
    "z3d": 1.1028045811246394,
    "element": "C",
    "radius": 18,
    "contacts": 15,
    "trajIndex": 11401
  },
  {
    "id": "L30",
    "x": 159.8586828377209,
    "y": -41.59108306646685,
    "x3d": 18.71999287181903,
    "y3d": -1.690949851400194,
    "z3d": 2.505648573170066,
    "element": "C",
    "radius": 18,
    "contacts": 0,
    "trajIndex": 11402
  },
  {
    "id": "L31",
    "x": 143.74465585765742,
    "y": -70.16565497027301,
    "x3d": 19.484558045711744,
    "y3d": 3.4646836254588123,
    "z3d": 1.4695531088339264,
    "element": "O",
    "radius": 18,
    "contacts": 233,
    "trajIndex": 11403
  },
  {
    "id": "L32",
    "x": 104.96847467050738,
    "y": 29.73258094036868,
    "x3d": 13.365085401369909,
    "y3d": -7.996179969210068,
    "z3d": -4.742913308726888,
    "element": "C",
    "radius": 18,
    "contacts": 330,
    "trajIndex": 11404
  },
  {
    "id": "L33",
    "x": 77.52337058690061,
    "y": 65.39441294378645,
    "x3d": 17.136695358591297,
    "y3d": -13.473473685615714,
    "z3d": -4.912135463025394,
    "element": "C",
    "radius": 18,
    "contacts": 404,
    "trajIndex": 11405
  },
  {
    "id": "L34",
    "x": 50.078266503293854,
    "y": 101.0562449472042,
    "x3d": 17.16512340818258,
    "y3d": -15.936830381656561,
    "z3d": -11.052803581047877,
    "element": "C",
    "radius": 18,
    "contacts": 403,
    "trajIndex": 11406
  },
  {
    "id": "L35",
    "x": 113.18520259031838,
    "y": 92.8395170273932,
    "x3d": 15.064769474407827,
    "y3d": -18.121645499461085,
    "z3d": -0.7066860737386177,
    "element": "C",
    "radius": 18,
    "contacts": 224,
    "trajIndex": 11407
  },
  {
    "id": "L36",
    "x": 204.46528733997695,
    "y": -35.65384172226425,
    "x3d": 22.963564412824443,
    "y3d": -5.07746198114968,
    "z3d": 4.896613869079921,
    "element": "N",
    "radius": 18,
    "contacts": 4,
    "trajIndex": 11408
  },
  {
    "id": "L37",
    "x": 231.91039142358375,
    "y": -71.315673725682,
    "x3d": 28.701956190357308,
    "y3d": -2.744195314837705,
    "z3d": 6.271889083668859,
    "element": "C",
    "radius": 18,
    "contacts": 0,
    "trajIndex": 11409
  },
  {
    "id": "L38",
    "x": 196.248559420166,
    "y": -98.76077780928875,
    "x3d": 32.77205888894196,
    "y3d": -7.906291343106034,
    "z3d": 6.70573448158829,
    "element": "C",
    "radius": 18,
    "contacts": 26,
    "trajIndex": 11410
  },
  {
    "id": "L39",
    "x": 202.18580076436857,
    "y": -143.36738231154484,
    "x3d": 38.32508543337435,
    "y3d": -6.618921759623161,
    "z3d": 5.599329535049854,
    "element": "O",
    "radius": 18,
    "contacts": 35,
    "trajIndex": 11411
  },
  {
    "id": "L40",
    "x": 267.57222342700146,
    "y": -43.87056964207526,
    "x3d": 28.619424322526278,
    "y3d": 0.6205171990565428,
    "z3d": 11.987158102325953,
    "element": "C",
    "radius": 18,
    "contacts": 18,
    "trajIndex": 11412
  },
  {
    "id": "L41",
    "x": 303.23405543041923,
    "y": -16.425465558468506,
    "x3d": 33.413782732817005,
    "y3d": 5.018351048410085,
    "z3d": 12.469951215174717,
    "element": "C",
    "radius": 18,
    "contacts": 46,
    "trajIndex": 11413
  },
  {
    "id": "L42",
    "x": 297.29681408621667,
    "y": 28.181138943787563,
    "x3d": 33.194909559817866,
    "y3d": 10.093479351451215,
    "z3d": 9.171196830078262,
    "element": "C",
    "radius": 18,
    "contacts": 111,
    "trajIndex": 11414
  },
  {
    "id": "L43",
    "x": 344.8331287780396,
    "y": -33.586965977117806,
    "x3d": 38.098925337454794,
    "y3d": 4.141946027938833,
    "z3d": 16.19477563799237,
    "element": "C",
    "radius": 18,
    "contacts": 128,
    "trajIndex": 11415
  },
  {
    "id": "L44",
    "x": 332.9586460896344,
    "y": 55.62624302739428,
    "x3d": 37.5656525838609,
    "y3d": 14.218810420567381,
    "z3d": 9.596258588026355,
    "element": "C",
    "radius": 18,
    "contacts": 299,
    "trajIndex": 11416
  },
  {
    "id": "L45",
    "x": 380.4949607814574,
    "y": -6.141861893511065,
    "x3d": 42.47197338287052,
    "y3d": 8.27390498453967,
    "z3d": 16.61759051653921,
    "element": "C",
    "radius": 18,
    "contacts": 314,
    "trajIndex": 11417
  },
  {
    "id": "L46",
    "x": 374.5577194372548,
    "y": 38.464742608745006,
    "x3d": 42.20098508811974,
    "y3d": 13.305593865082418,
    "z3d": 13.32121783670396,
    "element": "C",
    "radius": 18,
    "contacts": 286,
    "trajIndex": 11418
  },
  {
    "id": "L47",
    "x": 164.62396832792632,
    "y": -90.03992560043693,
    "x3d": 31.34386664805559,
    "y3d": -12.842810224304154,
    "z3d": 7.9344073275994225,
    "element": "O",
    "radius": 18,
    "contacts": 303,
    "trajIndex": 11419
  },
  {
    "id": "L48",
    "x": -319.87062537734363,
    "y": 85.04801076850605,
    "x3d": -47.22474573180085,
    "y3d": 5.720812402729925,
    "z3d": -3.9789320182755,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L49",
    "x": -328.08735329715466,
    "y": 21.941074681481574,
    "x3d": -50.0,
    "y3d": 10.240731125187606,
    "z3d": -9.513484216618714,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L50",
    "x": -256.76368929031923,
    "y": 76.83128284869514,
    "x3d": -46.89960091981335,
    "y3d": 3.4429314747446726,
    "z3d": -11.360248363135213,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L51",
    "x": -203.21231228922468,
    "y": 61.260661453493505,
    "x3d": -37.053296286836314,
    "y3d": 2.569059835685639,
    "z3d": -5.693741630033109,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L52",
    "x": -165.4836045406605,
    "y": -51.66207590096236,
    "x3d": -29.606511787651762,
    "y3d": 9.877752143048713,
    "z3d": -11.885818543583179,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L53",
    "x": -201.14543654407828,
    "y": -79.10717998456911,
    "x3d": -31.0688510228591,
    "y3d": 17.374995097354258,
    "z3d": -4.796742737479241,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L54",
    "x": -254.29417118722355,
    "y": -23.06174179728986,
    "x3d": -24.129823485372455,
    "y3d": 14.557325234265281,
    "z3d": -3.731879511014744,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L55",
    "x": -236.80726854749602,
    "y": -106.55228406817588,
    "x3d": -30.258266508455733,
    "y3d": 16.73964670747352,
    "z3d": 5.2105243178087495,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L56",
    "x": -295.83985296820657,
    "y": -38.84000610915433,
    "x3d": -34.70561791282817,
    "y3d": 10.806474782285907,
    "z3d": 3.4304454773221873,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L57",
    "x": -329.63879529373577,
    "y": -170.3872647342277,
    "x3d": -32.39771749359818,
    "y3d": 7.0949030341389285,
    "z3d": 12.551164315208203,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L58",
    "x": -252.37788994269755,
    "y": -160.10366106927026,
    "x3d": -25.673727061024614,
    "y3d": 3.9213205560251994,
    "z3d": 13.647836337149254,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L59",
    "x": -109.73056192902652,
    "y": -50.323244734843215,
    "x3d": -19.233583171244213,
    "y3d": 7.791705065613986,
    "z3d": -8.657978364468235,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L60",
    "x": -44.34413926639359,
    "y": 49.173567934626384,
    "x3d": -17.912378400900675,
    "y3d": -2.462221697870055,
    "z3d": -2.03882291254891,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L61",
    "x": -72.58386217575656,
    "y": 77.99435965805398,
    "x3d": -13.74879718369259,
    "y3d": -2.640078469980608,
    "z3d": -14.668912200009872,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L62",
    "x": -137.51726984642647,
    "y": 44.26810191788118,
    "x3d": -19.662817251521993,
    "y3d": -6.650470247375666,
    "z3d": -12.175674466736174,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L63",
    "x": -173.63752741220432,
    "y": 70.15872363161446,
    "x3d": -7.130246654915472,
    "y3d": -8.319627428349353,
    "z3d": -8.276324757192066,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L64",
    "x": -189.7863876042384,
    "y": 164.37579186469068,
    "x3d": -7.550938762595958,
    "y3d": -17.299778322067716,
    "z3d": -13.854990373523124,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L65",
    "x": -143.41581568830716,
    "y": 147.36605112552934,
    "x3d": -8.923043421348085,
    "y3d": -11.324543329768442,
    "z3d": -18.386157389221562,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L66",
    "x": -198.00311552404946,
    "y": 101.2688557776662,
    "x3d": -14.712495024512117,
    "y3d": -15.779428684638322,
    "z3d": -16.08138987423521,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L67",
    "x": -63.572515430189334,
    "y": 147.94233602506867,
    "x3d": -10.741309101363946,
    "y3d": -17.275393945038225,
    "z3d": -3.547756204248445,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L68",
    "x": -112.75267369515993,
    "y": 143.37361050732568,
    "x3d": -13.933634833316258,
    "y3d": -11.150317271322155,
    "z3d": -0.3898395851589673,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L69",
    "x": -70.8056151302137,
    "y": 91.65183058141857,
    "x3d": -17.92892034896943,
    "y3d": -15.412714885338843,
    "z3d": -5.324746165264542,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L70",
    "x": 9.207237734700819,
    "y": 33.60294653942487,
    "x3d": -7.839052857277075,
    "y3d": -3.3936630742649463,
    "z3d": 0.5881701403775305,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L71",
    "x": 46.93594548326499,
    "y": -79.31979081503107,
    "x3d": -0.07320048573018512,
    "y3d": 5.810927169309441,
    "z3d": -3.316793599681695,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L72",
    "x": -51.4297430859797,
    "y": -43.14201188855125,
    "x3d": -2.0946774507619454,
    "y3d": 4.0290602645528955,
    "z3d": 9.716855841855812,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L73",
    "x": 11.274113479847223,
    "y": -106.76489489863782,
    "x3d": -5.9805523813221555,
    "y3d": 9.413003799523297,
    "z3d": 5.805748616867159,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L74",
    "x": -81.55741326639246,
    "y": -170.59987556468965,
    "x3d": 6.2938970979694755,
    "y3d": 16.451578330389356,
    "z3d": 11.124517920121852,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L75",
    "x": 102.68898809489897,
    "y": -77.98095964891192,
    "x3d": 9.702057101727384,
    "y3d": 3.726905190465029,
    "z3d": 0.08365137534875443,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L76",
    "x": 168.0754107575319,
    "y": 21.515853020557664,
    "x3d": 11.945036962276815,
    "y3d": -7.7843450885054315,
    "z3d": 4.635658499402995,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L77",
    "x": 143.70982264568687,
    "y": 52.625985166609325,
    "x3d": 14.981008342255334,
    "y3d": -4.95668614451178,
    "z3d": -8.009683062752169,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L78",
    "x": 63.72540952234773,
    "y": 11.73231490942744,
    "x3d": 8.98780111671108,
    "y3d": -9.169206431627197,
    "z3d": -6.142138816270585,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L79",
    "x": 38.782022611721075,
    "y": 42.50100871754578,
    "x3d": 21.605718499042407,
    "y3d": -12.355849548043984,
    "z3d": -3.8445871071760576,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L80",
    "x": 22.633162419687093,
    "y": 136.718076950622,
    "x3d": 19.991180114324848,
    "y3d": -19.717077275797383,
    "z3d": -11.338356139430202,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L81",
    "x": 82.12871245889765,
    "y": 132.64372528430803,
    "x3d": 18.756240103850526,
    "y3d": -12.719983909294891,
    "z3d": -14.134713024972593,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L82",
    "x": 14.41643449987609,
    "y": 73.61114086359746,
    "x3d": 12.820426377964763,
    "y3d": -17.241452455492684,
    "z3d": -12.401615889173005,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L83",
    "x": 148.84703459373614,
    "y": 120.28462111099999,
    "x3d": 10.561547851308601,
    "y3d": -19.31943129680925,
    "z3d": -1.5304130100513866,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L84",
    "x": 90.29179836407769,
    "y": 131.58086500257272,
    "x3d": 15.414295693650292,
    "y3d": -16.69619661729994,
    "z3d": 3.7797233013669307,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L85",
    "x": 144.7726829274222,
    "y": 60.78907107178944,
    "x3d": 17.696230262847628,
    "y3d": -22.033564635045792,
    "z3d": -1.1209257507538173,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L86",
    "x": 216.4280475715982,
    "y": 7.726941669010969,
    "x3d": 22.152306317913805,
    "y3d": -9.249560171149094,
    "z3d": 5.985129086557938,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L87",
    "x": 259.3554955071905,
    "y": -106.97750572909979,
    "x3d": 30.016150065041813,
    "y3d": -0.06639900984206039,
    "z3d": 2.590668669524653,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L88",
    "x": 166.52396876095077,
    "y": -170.8124863951514,
    "x3d": 40.421280242382686,
    "y3d": -10.267040032652352,
    "z3d": 6.001690019877241,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L89",
    "x": 244.67881920076084,
    "y": -5.129221666895787,
    "x3d": 24.52249379042259,
    "y3d": 2.9551622232367967,
    "z3d": 12.48859535909544,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L90",
    "x": 295.0173275106083,
    "y": -79.53240164549302,
    "x3d": 28.81007284974444,
    "y3d": -2.3770111113072994,
    "z3d": 15.659849385250975,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L91",
    "x": 255.69774073859622,
    "y": 45.342639362436955,
    "x3d": 29.586310848966313,
    "y3d": 10.874045403048083,
    "z3d": 6.25816003260308,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L92",
    "x": 350.7703701222422,
    "y": -78.1935704793739,
    "x3d": 38.38275739768532,
    "y3d": 0.2483255650109937,
    "z3d": 18.802440759848913,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L93",
    "x": 327.02140474543177,
    "y": 100.23284752965036,
    "x3d": 37.346271599178976,
    "y3d": 18.142817588783206,
    "z3d": 7.026704455993325,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L94",
    "x": 422.09403412907767,
    "y": -23.30336231216031,
    "x3d": 46.09754169671032,
    "y3d": 7.5697451007916445,
    "z3d": 19.51630324404375,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  },
  {
    "id": "L95",
    "x": 410.2195514406726,
    "y": 65.90984669235178,
    "x3d": 45.60386446518117,
    "y3d": 16.516399512455937,
    "z3d": 13.653745821176408,
    "element": "H",
    "radius": 10,
    "contacts": 0,
    "trajIndex": -1
  }
];
        const CONTACTS = [
  {
    "id": "MET364",
    "frequency": 0.7195,
    "residueType": "MET",
    "ligandAtom": "L1",
    "residue": "MET364",
    "avgDistance": 0.33587609976530075,
    "pixelDistance": 321.38174655469544,
    "isTop3": true
  },
  {
    "id": "GLY174",
    "frequency": 0.616,
    "residueType": "GLY",
    "ligandAtom": "L26",
    "residue": "GLY174",
    "avgDistance": 0.3292710979779562,
    "pixelDistance": 319.6251834767662,
    "isTop3": true
  },
  {
    "id": "MET362",
    "frequency": 0.5366,
    "residueType": "MET",
    "ligandAtom": "L11",
    "residue": "MET362",
    "avgDistance": 0.3377023130655289,
    "pixelDistance": 321.85984094024514,
    "isTop3": true
  },
  {
    "id": "ARG365",
    "frequency": 0.48125,
    "residueType": "ARG",
    "ligandAtom": "L1",
    "residue": "ARG365",
    "avgDistance": 0.3327596900718553,
    "pixelDistance": 320.55840401890174,
    "isTop3": false
  },
  {
    "id": "PRO363",
    "frequency": 0.44283333333333336,
    "residueType": "PRO",
    "ligandAtom": "L12",
    "residue": "PRO363",
    "avgDistance": 0.33147019644578296,
    "pixelDistance": 320.21490131529674,
    "isTop3": false
  },
  {
    "id": "ARG152",
    "frequency": 0.43166666666666664,
    "residueType": "ARG",
    "ligandAtom": "L27",
    "residue": "ARG152",
    "avgDistance": 0.32155240376790367,
    "pixelDistance": 317.5145550329437,
    "isTop3": false
  },
  {
    "id": "HIS175",
    "frequency": 0.3785,
    "residueType": "HIS",
    "ligandAtom": "L8",
    "residue": "HIS175",
    "avgDistance": 0.3342202798773845,
    "pixelDistance": 320.9454755104363,
    "isTop3": false
  },
  {
    "id": "PRO242",
    "frequency": 0.312,
    "residueType": "PRO",
    "ligandAtom": "L44",
    "residue": "PRO242",
    "avgDistance": 0.3380148373544216,
    "pixelDistance": 321.94133999857667,
    "isTop3": false
  },
  {
    "id": "THR172",
    "frequency": 0.283,
    "residueType": "THR",
    "ligandAtom": "L34",
    "residue": "THR172",
    "avgDistance": 0.3361445963382721,
    "pixelDistance": 321.4522383322258,
    "isTop3": false
  },
  {
    "id": "VAL247",
    "frequency": 0.262,
    "residueType": "VAL",
    "ligandAtom": "L31",
    "residue": "VAL247",
    "avgDistance": 0.3350079581141472,
    "pixelDistance": 321.15334447009866,
    "isTop3": false
  },
  {
    "id": "PHE278",
    "frequency": 0.26,
    "residueType": "PHE",
    "ligandAtom": "L1",
    "residue": "PHE278",
    "avgDistance": 0.33397022634744644,
    "pixelDistance": 320.87935868804993,
    "isTop3": false
  },
  {
    "id": "VAL344",
    "frequency": 0.243,
    "residueType": "VAL",
    "ligandAtom": "L18",
    "residue": "VAL344",
    "avgDistance": 0.341672678788503,
    "pixelDistance": 322.88844527474066,
    "isTop3": false
  },
  {
    "id": "TYR323",
    "frequency": 0.161,
    "residueType": "TYR",
    "ligandAtom": "L10",
    "residue": "TYR323",
    "avgDistance": 0.32863684991995495,
    "pixelDistance": 319.45417101944497,
    "isTop3": false
  },
  {
    "id": "ARG176",
    "frequency": 0.141,
    "residueType": "ARG",
    "ligandAtom": "L34",
    "residue": "ARG176",
    "avgDistance": 0.34129246324300766,
    "pixelDistance": 322.79057199674133,
    "isTop3": false
  },
  {
    "id": "SER346",
    "frequency": 0.129,
    "residueType": "SER",
    "ligandAtom": "L35",
    "residue": "SER346",
    "avgDistance": 0.3362237215042114,
    "pixelDistance": 321.47299879260646,
    "isTop3": false
  }
];
        const LIGAND_BONDS = [
  [
    "L0",
    "L1"
  ],
  [
    "L0",
    "L2"
  ],
  [
    "L0",
    "L3"
  ],
  [
    "L3",
    "L4"
  ],
  [
    "L4",
    "L5"
  ],
  [
    "L4",
    "L7"
  ],
  [
    "L5",
    "L6"
  ],
  [
    "L5",
    "L12"
  ],
  [
    "L7",
    "L8"
  ],
  [
    "L8",
    "L9"
  ],
  [
    "L9",
    "L10"
  ],
  [
    "L9",
    "L11"
  ],
  [
    "L12",
    "L13"
  ],
  [
    "L13",
    "L14"
  ],
  [
    "L13",
    "L16"
  ],
  [
    "L14",
    "L15"
  ],
  [
    "L14",
    "L20"
  ],
  [
    "L16",
    "L17"
  ],
  [
    "L17",
    "L18"
  ],
  [
    "L17",
    "L19"
  ],
  [
    "L20",
    "L21"
  ],
  [
    "L21",
    "L22"
  ],
  [
    "L21",
    "L24"
  ],
  [
    "L22",
    "L23"
  ],
  [
    "L22",
    "L28"
  ],
  [
    "L24",
    "L25"
  ],
  [
    "L25",
    "L26"
  ],
  [
    "L25",
    "L27"
  ],
  [
    "L28",
    "L29"
  ],
  [
    "L29",
    "L30"
  ],
  [
    "L29",
    "L32"
  ],
  [
    "L30",
    "L31"
  ],
  [
    "L30",
    "L36"
  ],
  [
    "L32",
    "L33"
  ],
  [
    "L33",
    "L34"
  ],
  [
    "L33",
    "L35"
  ],
  [
    "L36",
    "L37"
  ],
  [
    "L37",
    "L38"
  ],
  [
    "L37",
    "L40"
  ],
  [
    "L38",
    "L39"
  ],
  [
    "L38",
    "L47"
  ],
  [
    "L40",
    "L41"
  ],
  [
    "L41",
    "L42"
  ],
  [
    "L41",
    "L43"
  ],
  [
    "L42",
    "L44"
  ],
  [
    "L43",
    "L45"
  ],
  [
    "L44",
    "L46"
  ],
  [
    "L45",
    "L46"
  ],
  [
    "L2",
    "L48"
  ],
  [
    "L2",
    "L49"
  ],
  [
    "L2",
    "L50"
  ],
  [
    "L3",
    "L51"
  ],
  [
    "L4",
    "L52"
  ],
  [
    "L7",
    "L53"
  ],
  [
    "L7",
    "L54"
  ],
  [
    "L8",
    "L55"
  ],
  [
    "L8",
    "L56"
  ],
  [
    "L11",
    "L57"
  ],
  [
    "L11",
    "L58"
  ],
  [
    "L12",
    "L59"
  ],
  [
    "L13",
    "L60"
  ],
  [
    "L16",
    "L61"
  ],
  [
    "L16",
    "L62"
  ],
  [
    "L17",
    "L63"
  ],
  [
    "L18",
    "L64"
  ],
  [
    "L18",
    "L65"
  ],
  [
    "L18",
    "L66"
  ],
  [
    "L19",
    "L67"
  ],
  [
    "L19",
    "L68"
  ],
  [
    "L19",
    "L69"
  ],
  [
    "L20",
    "L70"
  ],
  [
    "L21",
    "L71"
  ],
  [
    "L24",
    "L72"
  ],
  [
    "L24",
    "L73"
  ],
  [
    "L27",
    "L74"
  ],
  [
    "L28",
    "L75"
  ],
  [
    "L29",
    "L76"
  ],
  [
    "L32",
    "L77"
  ],
  [
    "L32",
    "L78"
  ],
  [
    "L33",
    "L79"
  ],
  [
    "L34",
    "L80"
  ],
  [
    "L34",
    "L81"
  ],
  [
    "L34",
    "L82"
  ],
  [
    "L35",
    "L83"
  ],
  [
    "L35",
    "L84"
  ],
  [
    "L35",
    "L85"
  ],
  [
    "L36",
    "L86"
  ],
  [
    "L37",
    "L87"
  ],
  [
    "L39",
    "L88"
  ],
  [
    "L40",
    "L89"
  ],
  [
    "L40",
    "L90"
  ],
  [
    "L42",
    "L91"
  ],
  [
    "L43",
    "L92"
  ],
  [
    "L44",
    "L93"
  ],
  [
    "L45",
    "L94"
  ],
  [
    "L46",
    "L95"
  ]
];
        class ContactMap {
            constructor() {
                try {
                    this.canvas = document.getElementById('canvas');
                    this.ctx = this.canvas.getContext('2d');
                    this.tooltip = document.getElementById('tooltip');
                    this.showConnections = true;
                    this.showHydrogens = true;  
                    this.distanceLocked = true;
                    this.is3DMode = false;
                    this.isDragging = false;
                    this.isPanning = false;
                    this.dragTarget = null;
                    this.wheelRotation = 0;
                    this.offset = {x: 0, y: 0};
                    this.panStart = {x: 0, y: 0};
                    this.viewOffset = {x: 0, y: 0};
                    this.zoom = 1.0;
                    this.centerX = 600;
                    this.centerY = 400;
                    
                    // 3D rotation angles
                    this.rotationX = 0;
                    this.rotationY = 0;
                    this.rotationZ = 0;
                    
                    // Initialize data
                    this.initializeData();
                    
                    // Store initial state for reset
                    this.saveInitialState();
                    
                    this.initEvents();
                    this.animate();
                } catch(e) {
                    console.error('Error initializing ContactMap:', e);
                }
            }
            
            initializeData() {
                // Calculate ligand center
                let ligandCenterX = 0;
                let ligandCenterY = 0;
                if (LIGAND_ATOMS.length > 0) {
                    LIGAND_ATOMS.forEach(atom => {
                        ligandCenterX += atom.x;
                        ligandCenterY += atom.y;
                    });
                    ligandCenterX /= LIGAND_ATOMS.length;
                    ligandCenterY /= LIGAND_ATOMS.length;
                }
                
                // Initialize ligand atoms with proper 3D coordinates
                this.ligandAtoms = LIGAND_ATOMS.map(atom => {
                    // Use provided 3D coordinates or generate from 2D
                    let x3d = atom.x3d;
                    let y3d = atom.y3d;
                    let z3d = atom.z3d;
                    
                    // If 3D coordinates are missing or invalid, use 2D coordinates
                    if (x3d === undefined || y3d === undefined || z3d === undefined) {
                        x3d = atom.x;
                        y3d = atom.y;
                        z3d = 0;
                    }
                    
                    return {
                        ...atom,
                        x: atom.x + this.centerX,
                        y: atom.y + this.centerY,
                        x3d: x3d * 3,  // Scale up 3D coordinates for better visibility
                        y3d: y3d * 3,
                        z3d: z3d * 3,
                        fixed: true
                    };
                });
                
                // Group contacts by ligand atom to handle overlaps
                const contactsByLigandAtom = {};
                CONTACTS.forEach(contact => {
                    const atomId = contact.ligandAtom || 'L0';
                    if (!contactsByLigandAtom[atomId]) {
                        contactsByLigandAtom[atomId] = [];
                    }
                    contactsByLigandAtom[atomId].push(contact);
                });
                
                // Initialize contacts with proper alignment and spacing
                this.contacts = [];
                
                // Process each group of contacts
                Object.keys(contactsByLigandAtom).forEach(ligandAtomId => {
                    const contactGroup = contactsByLigandAtom[ligandAtomId];
                    const ligandAtom = this.ligandAtoms.find(a => a.id === ligandAtomId);
                    
                    if (ligandAtom) {
                        // Calculate vector from ligand center to this atom
                        const atomRelX = ligandAtom.x - this.centerX;
                        const atomRelY = ligandAtom.y - this.centerY;
                        
                        // Calculate base angle from center through atom
                        let baseAngle = Math.atan2(atomRelY, atomRelX);
                        const baseDirLength = Math.sqrt(atomRelX * atomRelX + atomRelY * atomRelY);
                        
                        // If atom is at center, distribute evenly
                        if (baseDirLength < 1) {
                            baseAngle = 0;
                        }
                        
                        // For multiple residues on same atom, create larger angular offset for TOP3
                        const spreadAngle = Math.PI / 6; // 30 degrees total spread
                        
                        contactGroup.forEach((contact, groupIdx) => {
                            // Calculate offset angle for multiple contacts
                            let angle = baseAngle;
                            if (contactGroup.length > 1) {
                                // Give TOP3 contacts more space
                                const spreadFactor = contact.isTop3 ? 1.5 : 1.0;
                                const offset = (groupIdx - (contactGroup.length - 1) / 2) * 
                                             (spreadAngle * spreadFactor / Math.max(1, contactGroup.length - 1));
                                angle = baseAngle + offset;
                            }
                            
                            // Adjust distance based on whether it's TOP3
                            const baseDistance = contact.pixelDistance || 250;
                            const distance = contact.isTop3 ? baseDistance * 0.7 : baseDistance * 0.5;
                            
                            // Calculate position along the extension line
                            const dirX = Math.cos(angle);
                            const dirY = Math.sin(angle);
                            
                            // Position residue along the line from ligand center through ligand atom
                            const x = ligandAtom.x + dirX * distance;
                            const y = ligandAtom.y + dirY * distance;
                            
                            this.contacts.push({
                                ...contact,
                                x: x,
                                y: y,
                                angle: angle,
                                radius: distance,
                                fixed: false,
                                initialX: x,
                                initialY: y,
                                pixelDistance: distance
                            });
                        });
                    } else {
                        // Fallback positioning
                        contactGroup.forEach((contact, idx) => {
                            const angle = (this.contacts.length * 2 * Math.PI / CONTACTS.length) - Math.PI/2;
                            const baseDistance = contact.pixelDistance || 250;
                            const distance = contact.isTop3 ? baseDistance * 0.7 : baseDistance * 0.5;
                            const x = this.centerX + Math.cos(angle) * distance;
                            const y = this.centerY + Math.sin(angle) * distance;
                            
                            this.contacts.push({
                                ...contact,
                                x: x,
                                y: y,
                                angle: angle,
                                radius: distance,
                                fixed: false,
                                initialX: x,
                                initialY: y,
                                pixelDistance: distance
                            });
                        });
                    }
                });
                
                this.ligandBonds = LIGAND_BONDS;
            }
            
            saveInitialState() {
                this.initialState = {
                    ligandAtoms: this.ligandAtoms.map(atom => ({...atom})),
                    contacts: this.contacts.map(contact => ({...contact})),
                    zoom: 1.0,
                    viewOffset: {x: 0, y: 0},
                    rotationX: 0,
                    rotationY: 0,
                    rotationZ: 0
                };
            }
            
            restoreInitialState() {
                if (!this.initialState) return;
                
                this.ligandAtoms = this.initialState.ligandAtoms.map(atom => ({...atom}));
                this.contacts = this.initialState.contacts.map(contact => ({...contact}));
                this.zoom = this.initialState.zoom;
                this.viewOffset = {...this.initialState.viewOffset};
                this.rotationX = this.initialState.rotationX;
                this.rotationY = this.initialState.rotationY;
                this.rotationZ = this.initialState.rotationZ;
                this.wheelRotation = 0;
            }
            
            initEvents() {
                this.canvas.addEventListener('mousedown', e => this.onMouseDown(e));
                this.canvas.addEventListener('mousemove', e => this.onMouseMove(e));
                this.canvas.addEventListener('mouseup', e => this.onMouseUp(e));
                this.canvas.addEventListener('mouseleave', e => this.onMouseLeave(e));
                this.canvas.addEventListener('wheel', e => this.onWheel(e));
                this.canvas.addEventListener('contextmenu', e => e.preventDefault());
            }
            
            getMousePos(e) {
                const rect = this.canvas.getBoundingClientRect();
                const x = (e.clientX - rect.left - this.viewOffset.x) / this.zoom;
                const y = (e.clientY - rect.top - this.viewOffset.y) / this.zoom;
                return { x, y };
            }
            
            onMouseDown(e) {
                const pos = this.getMousePos(e);
                
                if (e.button === 1 || (e.button === 0 && e.shiftKey)) {
                    // Middle click or Shift+click for panning
                    this.isPanning = true;
                    this.panStart = { x: e.clientX - this.viewOffset.x, y: e.clientY - this.viewOffset.y };
                    this.canvas.style.cursor = 'move';
                    e.preventDefault();
                } else if (e.button === 0) {
                    if (this.is3DMode) {
                        // 3D mode rotation
                        this.isDragging = true;
                        this.dragStart = { x: e.clientX, y: e.clientY };
                        this.startRotationX = this.rotationX;
                        this.startRotationY = this.rotationY;
                        this.canvas.style.cursor = 'grabbing';
                    } else {
                        // 2D mode
                        this.dragTarget = this.findElementAt(pos.x, pos.y);
                        
                        if (this.dragTarget && !this.dragTarget.fixed && this.dragTarget.radius) {
                            this.isDragging = true;
                            const ligandAtom = this.ligandAtoms.find(a => a.id === this.dragTarget.ligandAtom);
                            this.rotationCenter = ligandAtom ? { x: ligandAtom.x, y: ligandAtom.y } : { x: this.centerX, y: this.centerY };
                            
                            if (this.distanceLocked) {
                                const dx = pos.x - this.rotationCenter.x;
                                const dy = pos.y - this.rotationCenter.y;
                                this.dragStartAngle = Math.atan2(dy, dx);
                                this.targetStartAngle = Math.atan2(this.dragTarget.y - this.rotationCenter.y, 
                                                                 this.dragTarget.x - this.rotationCenter.x);
                            } else {
                                this.offset = { x: pos.x - this.dragTarget.x, y: pos.y - this.dragTarget.y };
                            }
                            this.canvas.style.cursor = 'grabbing';
                        } else if (!this.dragTarget) {
                            this.isPanning = true;
                            this.panStart = { x: e.clientX - this.viewOffset.x, y: e.clientY - this.viewOffset.y };
                            this.canvas.style.cursor = 'move';
                        }
                    }
                }
            }
            
            onMouseMove(e) {
                if (this.isPanning) {
                    this.viewOffset.x = e.clientX - this.panStart.x;
                    this.viewOffset.y = e.clientY - this.panStart.y;
                } else if (this.isDragging) {
                    if (this.is3DMode) {
                        const dx = e.clientX - this.dragStart.x;
                        const dy = e.clientY - this.dragStart.y;
                        this.rotationY = this.startRotationY + dx * 0.01;
                        this.rotationX = this.startRotationX + dy * 0.01;
                        this.rotationX = Math.max(-Math.PI, Math.min(Math.PI, this.rotationX));
                    } else if (this.dragTarget) {
                        const pos = this.getMousePos(e);
                        
                        if (this.distanceLocked) {
                            const dx = pos.x - this.rotationCenter.x;
                            const dy = pos.y - this.rotationCenter.y;
                            const currentAngle = Math.atan2(dy, dx);
                            const angleDiff = currentAngle - this.dragStartAngle;
                            
                            const dist = Math.sqrt(Math.pow(this.dragTarget.x - this.rotationCenter.x, 2) + 
                                                 Math.pow(this.dragTarget.y - this.rotationCenter.y, 2));
                            
                            const newAngle = this.targetStartAngle + angleDiff;
                            this.dragTarget.x = this.rotationCenter.x + dist * Math.cos(newAngle);
                            this.dragTarget.y = this.rotationCenter.y + dist * Math.sin(newAngle);
                            this.dragTarget.angle = newAngle;
                        } else {
                            this.dragTarget.x = pos.x - this.offset.x;
                            this.dragTarget.y = pos.y - this.offset.y;
                            
                            const dx = this.dragTarget.x - this.rotationCenter.x;
                            const dy = this.dragTarget.y - this.rotationCenter.y;
                            this.dragTarget.angle = Math.atan2(dy, dx);
                            this.dragTarget.radius = Math.sqrt(dx * dx + dy * dy);
                        }
                    }
                } else {
                    const pos = this.getMousePos(e);
                    this.handleHover(pos.x, pos.y, e);
                }
            }
            
            onMouseUp(e) {
                this.isDragging = false;
                this.isPanning = false;
                this.dragTarget = null;
                this.canvas.style.cursor = 'grab';
            }
            
            onMouseLeave(e) {
                this.onMouseUp(e);
                this.hideTooltip();
            }
            
            onWheel(e) {
                e.preventDefault();
                const delta = e.deltaY > 0 ? 0.9 : 1.1;
                const newZoom = this.zoom * delta;
                
                if (newZoom >= 0.3 && newZoom <= 3) {
                    const rect = this.canvas.getBoundingClientRect();
                    const mouseX = e.clientX - rect.left;
                    const mouseY = e.clientY - rect.top;
                    
                    const worldX = (mouseX - this.viewOffset.x) / this.zoom;
                    const worldY = (mouseY - this.viewOffset.y) / this.zoom;
                    
                    this.zoom = newZoom;
                    
                    this.viewOffset.x = mouseX - worldX * this.zoom;
                    this.viewOffset.y = mouseY - worldY * this.zoom;
                }
            }
            
            findElementAt(x, y) {
                // Check contacts first (smaller hit areas for TOP3)
                for (let contact of this.contacts) {
                    const hitRadius = contact.isTop3 ? 50 : 40;
                    const dist = Math.sqrt((x - contact.x) ** 2 + (y - contact.y) ** 2);
                    if (dist < hitRadius) return contact;
                }
                for (let atom of this.ligandAtoms) {
                    const dist = Math.sqrt((x - atom.x) ** 2 + (y - atom.y) ** 2);
                    if (dist < atom.radius) return atom;
                }
                return null;
            }
            
            handleHover(x, y, e) {
                const element = this.findElementAt(x, y);
                if (element) {
                    this.showTooltip(element, e.clientX, e.clientY);
                    this.canvas.style.cursor = element.fixed ? 'default' : 'grab';
                } else {
                    this.hideTooltip();
                    this.canvas.style.cursor = 'grab';
                }
            }
            
            showTooltip(element, x, y) {
                let content = '';
                if (element.frequency !== undefined) {
                    const distanceText = element.avgDistance ? ` | Distance: ${element.avgDistance.toFixed(2)}nm` : '';
                    const top3Text = element.isTop3 ? ' <span style="color: gold;">★ TOP 3</span>' : '';
                    content = `${element.residue}${top3Text}<br>Frequency: ${(element.frequency * 100).toFixed(1)}%${distanceText}`;
                } else {
                    content = `Ligand Atom: ${element.element}${element.id.substring(1)}<br>Contacts: ${element.contacts || 0}`;
                }
                this.tooltip.innerHTML = content;
                this.tooltip.style.display = 'block';
                this.tooltip.style.left = (x + 10) + 'px';
                this.tooltip.style.top = (y - 10) + 'px';
            }
            
            hideTooltip() {
                this.tooltip.style.display = 'none';
            }
            
            rotate3D(x, y, z, rotX, rotY, rotZ) {
                let x1 = x, y1 = y, z1 = z;
                const cosX = Math.cos(rotX);
                const sinX = Math.sin(rotX);
                const y2 = y1 * cosX - z1 * sinX;
                const z2 = y1 * sinX + z1 * cosX;
                const cosY = Math.cos(rotY);
                const sinY = Math.sin(rotY);
                const x3 = x1 * cosY + z2 * sinY;
                const z3 = -x1 * sinY + z2 * cosY;
                const cosZ = Math.cos(rotZ);
                const sinZ = Math.sin(rotZ);
                const x4 = x3 * cosZ - y2 * sinZ;
                const y4 = x3 * sinZ + y2 * cosZ;
                return { x: x4, y: y4, z: z3 };
            }
            
            project3D(x, y, z) {
                const focalLength = 800;
                const scale = focalLength / (focalLength + z);
                return {
                    x: x * scale,
                    y: y * scale,
                    scale: scale,
                    depth: z
                };
            }
            
            draw3DScene() {
                const elements3D = [];
                this.ligandAtoms.forEach(atom => {
                    const rotated = this.rotate3D(
                        atom.x3d * 2, atom.y3d * 2, atom.z3d * 2,
                        this.rotationX, this.rotationY, this.rotationZ
                    );
                    const projected = this.project3D(rotated.x, rotated.y, rotated.z);
                    elements3D.push({
                        type: 'atom',
                        data: atom,
                        x: this.centerX + projected.x,
                        y: this.centerY + projected.y,
                        z: projected.depth,
                        scale: projected.scale
                    });
                });
                
                this.ligandBonds.forEach(([id1, id2]) => {
                    const atom1 = this.ligandAtoms.find(a => a.id === id1);
                    const atom2 = this.ligandAtoms.find(a => a.id === id2);
                    if (atom1 && atom2) {
                        const rot1 = this.rotate3D(
                            atom1.x3d * 2, atom1.y3d * 2, atom1.z3d * 2,
                            this.rotationX, this.rotationY, this.rotationZ
                        );
                        const rot2 = this.rotate3D(
                            atom2.x3d * 2, atom2.y3d * 2, atom2.z3d * 2,
                            this.rotationX, this.rotationY, this.rotationZ
                        );
                        const proj1 = this.project3D(rot1.x, rot1.y, rot1.z);
                        const proj2 = this.project3D(rot2.x, rot2.y, rot2.z);
                        elements3D.push({
                            type: 'bond',
                            x1: this.centerX + proj1.x,
                            y1: this.centerY + proj1.y,
                            x2: this.centerX + proj2.x,
                            y2: this.centerY + proj2.y,
                            z: (proj1.depth + proj2.depth) / 2,
                            scale: (proj1.scale + proj2.scale) / 2,
                            atom1: atom1,
                            atom2: atom2
                        });
                    }
                });
                
                elements3D.sort((a, b) => b.z - a.z);
                
                elements3D.forEach(elem => {
                    if (elem.type === 'bond') {
                        this.draw3DBond(elem);
                    } else if (elem.type === 'atom') {
                        this.draw3DAtom(elem);
                    }
                });
                
                this.draw3DAxisIndicator();
            }
            
            draw3DBond(bond) {
                const width = 3 * bond.scale;
                const opacity = 0.3 + bond.scale * 0.5;
                this.ctx.globalAlpha = opacity;
                this.ctx.strokeStyle = '#34495e';
                this.ctx.lineWidth = width;
                if (bond.atom1.element === 'H' || bond.atom2.element === 'H') {
                    this.ctx.lineWidth = width * 0.5;
                }
                this.ctx.beginPath();
                this.ctx.moveTo(bond.x1, bond.y1);
                this.ctx.lineTo(bond.x2, bond.y2);
                this.ctx.stroke();
                this.ctx.globalAlpha = 1.0;
            }
            
            draw3DAtom(elem) {
                const atom = elem.data;
                const radius = atom.radius * elem.scale;
                const opacity = 0.7 + elem.scale * 0.3;
                
                this.ctx.globalAlpha = opacity;
                this.ctx.fillStyle = this.getAtomColor(atom.element);
                this.ctx.beginPath();
                this.ctx.arc(elem.x, elem.y, radius, 0, 2 * Math.PI);
                this.ctx.fill();
                
                const gradient = this.ctx.createRadialGradient(
                    elem.x - radius * 0.3, elem.y - radius * 0.3, 0,
                    elem.x, elem.y, radius
                );
                gradient.addColorStop(0, 'rgba(255, 255, 255, 0.5)');
                gradient.addColorStop(0.5, 'rgba(255, 255, 255, 0.1)');
                gradient.addColorStop(1, 'rgba(0, 0, 0, 0.4)');
                this.ctx.fillStyle = gradient;
                this.ctx.fill();
                
                this.ctx.strokeStyle = elem.z > 0 ? 'rgba(255,255,255,0.8)' : 'rgba(0,0,0,0.3)';
                this.ctx.lineWidth = atom.element === 'H' ? 1 : 2;
                this.ctx.stroke();
                
                if (atom.element !== 'H' || elem.scale > 0.7) {
                    this.ctx.globalAlpha = opacity;
                    this.ctx.fillStyle = atom.element === 'H' ? '#333333' : '#ffffff';
                    this.ctx.font = `bold ${Math.floor(12 * elem.scale)}px Arial`;
                    this.ctx.textAlign = 'center';
                    this.ctx.textBaseline = 'middle';
                    const label = atom.element + atom.id.substring(1);
                    this.ctx.fillText(label, elem.x, elem.y);
                }
                this.ctx.globalAlpha = 1.0;
            }
            
            draw3DAxisIndicator() {
                const x = 80;
                const y = 80;
                const length = 40;
                const axes = [
                    { x: 1, y: 0, z: 0, color: '#e74c3c', label: 'X' },
                    { x: 0, y: 1, z: 0, color: '#2ecc71', label: 'Y' },
                    { x: 0, y: 0, z: 1, color: '#3498db', label: 'Z' }
                ];
                
                axes.forEach(axis => {
                    const rotated = this.rotate3D(
                        axis.x * length, axis.y * length, axis.z * length,
                        this.rotationX, this.rotationY, this.rotationZ
                    );
                    this.ctx.strokeStyle = axis.color;
                    this.ctx.lineWidth = 2;
                    this.ctx.beginPath();
                    this.ctx.moveTo(x, y);
                    this.ctx.lineTo(x + rotated.x, y + rotated.y);
                    this.ctx.stroke();
                    this.ctx.fillStyle = axis.color;
                    this.ctx.font = 'bold 12px Arial';
                    this.ctx.textAlign = 'center';
                    this.ctx.textBaseline = 'middle';
                    this.ctx.fillText(axis.label, x + rotated.x * 1.2, y + rotated.y * 1.2);
                });
                
                this.ctx.fillStyle = '#34495e';
                this.ctx.beginPath();
                this.ctx.arc(x, y, 3, 0, 2 * Math.PI);
                this.ctx.fill();
            }
            
            getAtomColor(element) {
                const colors = {
                    'C': '#27ae60',
                    'N': '#2196f3',
                    'O': '#e91e63',
                    'S': '#ffc107',
                    'H': '#ffffff',
                    'P': '#ff5722',
                    'F': '#9c27b0',
                    'Cl': '#4caf50',
                    'Br': '#795548',
                    'I': '#607d8b'
                };
                return colors[element] || '#7f8c8d';
            }
            
            getFrequencyColor(freq) {
                if (freq <= 0.3) {
                    const t = freq / 0.3;
                    const r = Math.floor(149 + (52 - 149) * t);
                    const g = Math.floor(165 + (152 - 165) * t);
                    const b = Math.floor(166 + (219 - 166) * t);
                    return `rgb(${r}, ${g}, ${b})`;
                } else if (freq <= 0.6) {
                    const t = (freq - 0.3) / 0.3;
                    const r = Math.floor(52 + (142 - 52) * t);
                    const g = Math.floor(152 + (68 - 152) * t);
                    const b = Math.floor(219 + (173 - 219) * t);
                    return `rgb(${r}, ${g}, ${b})`;
                } else {
                    const t = (freq - 0.6) / 0.4;
                    const r = Math.floor(142 + (231 - 142) * t);
                    const g = Math.floor(68 + (76 - 68) * t);
                    const b = Math.floor(173 + (60 - 173) * t);
                    return `rgb(${r}, ${g}, ${b})`;
                }
            }
            
            drawGrid() {
                this.ctx.strokeStyle = 'rgba(236, 240, 241, 0.3)';
                this.ctx.lineWidth = 1;
                const step = 50;
                const startX = -this.viewOffset.x / this.zoom;
                const startY = -this.viewOffset.y / this.zoom;
                const endX = (this.canvas.width - this.viewOffset.x) / this.zoom;
                const endY = (this.canvas.height - this.viewOffset.y) / this.zoom;
                
                for (let x = Math.floor(startX / step) * step; x <= endX; x += step) {
                    this.ctx.beginPath();
                    this.ctx.moveTo(x, startY);
                    this.ctx.lineTo(x, endY);
                    this.ctx.stroke();
                }
                for (let y = Math.floor(startY / step) * step; y <= endY; y += step) {
                    this.ctx.beginPath();
                    this.ctx.moveTo(startX, y);
                    this.ctx.lineTo(endX, y);
                    this.ctx.stroke();
                }
            }
            
            drawConnections() {
                this.contacts.forEach((contact, index) => {
                    const ligandAtom = this.ligandAtoms.find(a => a.id === contact.ligandAtom);
                    if (!ligandAtom) return;
                    
                    const color = this.getFrequencyColor(contact.frequency);
                    
                    let targetX = contact.x;
                    let targetY = contact.y;
                    
                    // Adjust target for TOP3 contacts
                    if (contact.isTop3) {
                        if (contact.contactAtomX !== undefined && contact.contactAtomY !== undefined) {
                            targetX = contact.contactAtomX;
                            targetY = contact.contactAtomY;
                        } else {
                            const angle = contact.angle || Math.atan2(ligandAtom.y - contact.y, ligandAtom.x - contact.x);
                            targetX = contact.x + Math.cos(angle) * 25;  // Reduced from 30
                            targetY = contact.y + Math.sin(angle) * 25;
                        }
                    }
                    
                    // Draw enhanced connections for TOP3
                    if (contact.isTop3) {
                        // Draw glow effect
                        this.ctx.save();
                        this.ctx.shadowBlur = 15;
                        this.ctx.shadowColor = color;
                        this.ctx.strokeStyle = color;
                        this.ctx.lineWidth = 4 + contact.frequency * 4;
                        this.ctx.globalAlpha = 0.3;
                        this.ctx.beginPath();
                        this.ctx.moveTo(ligandAtom.x, ligandAtom.y);
                        this.ctx.lineTo(targetX, targetY);
                        this.ctx.stroke();
                        this.ctx.restore();
                    }
                    
                    // Draw main connection
                    this.ctx.strokeStyle = color;
                    this.ctx.lineWidth = contact.isTop3 ? (3 + contact.frequency * 3) : (1.5 + contact.frequency * 2);
                    this.ctx.globalAlpha = contact.isTop3 ? 0.8 : (0.5 + contact.frequency * 0.3);
                    this.ctx.beginPath();
                    this.ctx.moveTo(ligandAtom.x, ligandAtom.y);
                    this.ctx.lineTo(targetX, targetY);
                    this.ctx.stroke();
                    this.ctx.globalAlpha = 1.0;
                });
            }
            
            drawLigand() {
                this.ctx.strokeStyle = '#34495e';
                this.ctx.lineWidth = 3;
                
                // 绘制键
                this.ligandBonds.forEach(([id1, id2]) => {
                    const atom1 = this.ligandAtoms.find(a => a.id === id1);
                    const atom2 = this.ligandAtoms.find(a => a.id === id2);
                    if (atom1 && atom2) {
                        // 如果隐藏氢原子模式开启，跳过涉及氢的键
                        if (!this.showHydrogens && (atom1.element === 'H' || atom2.element === 'H')) {
                            return;
                        }
                        
                        if (atom1.element === 'H' || atom2.element === 'H') {
                            this.ctx.lineWidth = 1.5;
                        } else {
                            this.ctx.lineWidth = 3;
                        }
                        this.ctx.beginPath();
                        this.ctx.moveTo(atom1.x, atom1.y);
                        this.ctx.lineTo(atom2.x, atom2.y);
                        this.ctx.stroke();
                    }
                });
                
                // 绘制原子
                this.ligandAtoms.forEach(atom => {
                    // 如果隐藏氢原子模式开启，跳过氢原子
                    if (!this.showHydrogens && atom.element === 'H') {
                        return;
                    }
                    
                    const color = this.getAtomColor(atom.element);
                    
                    this.ctx.fillStyle = color;
                    this.ctx.beginPath();
                    this.ctx.arc(atom.x, atom.y, atom.radius, 0, 2 * Math.PI);
                    this.ctx.fill();
                    
                    this.ctx.strokeStyle = atom.element === 'H' ? '#cccccc' : '#ffffff';
                    this.ctx.lineWidth = atom.element === 'H' ? 1 : 2;
                    this.ctx.stroke();
                    
                    const atomNumber = atom.id.substring(1);
                    let label = atom.element;
                    
                    if (atom.element === 'H') {
                        label = 'H' + atomNumber;
                        this.ctx.fillStyle = '#333333';
                        this.ctx.font = 'bold 9px Arial';
                    } else {
                        label = atom.element + atomNumber;
                        this.ctx.fillStyle = '#ffffff';
                        this.ctx.font = 'bold 12px Arial';
                    }
                    
                    this.ctx.textAlign = 'center';
                    this.ctx.textBaseline = 'middle';
                    this.ctx.fillText(label, atom.x, atom.y);
                });
            }
            
            drawContacts() {
                this.contacts.forEach((contact, i) => {
                    const color = this.getFrequencyColor(contact.frequency);
                    this.drawOrientedStructure(contact, color, i + 1);
                });
            }
            
            drawOrientedStructure(contact, color, index) {
                const x = contact.x, y = contact.y;
                const residueType = contact.residueType;
                
                const ligandAtom = this.ligandAtoms.find(a => a.id === contact.ligandAtom);
                if (!ligandAtom) {
                    this.drawSimpleLabel(x, y, residueType, contact.residue, color, index);
                    return;
                }
                
                const dx = ligandAtom.x - x;
                const dy = ligandAtom.y - y;
                const angle = Math.atan2(dy, dx);
                
                // Use contact.isTop3 instead of index
                if (contact.isTop3) {
                    this.drawOrientedAAStructure(x, y, residueType, contact.residue, color, index, angle, contact.frequency);
                } else {
                    this.drawSimpleLabel(x, y, residueType, contact.residue, color, index);
                }
            }
            
            drawOrientedAAStructure(x, y, residueType, residueName, color, index, angle, frequency) {
                // Calculate contact atom position (closer to structure)
                const contactAtomPos = {
                    x: x + Math.cos(angle) * 25,  // Reduced from 30
                    y: y + Math.sin(angle) * 25
                };
                
                const contactIdx = this.contacts.findIndex(c => c.residue === residueName);
                if (contactIdx >= 0) {
                    this.contacts[contactIdx].contactAtomX = contactAtomPos.x;
                    this.contacts[contactIdx].contactAtomY = contactAtomPos.y;
                }
                
                this.ctx.save();
                this.ctx.translate(x, y);
                this.ctx.rotate(angle);
                
                // Smaller backbone structure
                const backbone = [
                    {name: 'N', x: 22, y: 0, color: '#4169E1', radius: 11},    // Reduced sizes
                    {name: 'Ca', x: 7, y: 0, color: '#808080', radius: 11},
                    {name: 'C', x: -7, y: 0, color: '#808080', radius: 11},
                    {name: 'O1', x: -22, y: -7, color: '#DC143C', radius: 9},
                    {name: 'O2', x: -22, y: 7, color: '#DC143C', radius: 9}
                ];
                
                // Draw bonds with thinner lines
                this.ctx.strokeStyle = '#000000';
                this.ctx.lineWidth = 4;  // Reduced from 5
                
                this.ctx.beginPath();
                this.ctx.moveTo(backbone[0].x - backbone[0].radius * 0.8, backbone[0].y);
                this.ctx.lineTo(backbone[1].x + backbone[1].radius * 0.8, backbone[1].y);
                this.ctx.stroke();
                
                this.ctx.beginPath();
                this.ctx.moveTo(backbone[1].x - backbone[1].radius * 0.8, backbone[1].y);
                this.ctx.lineTo(backbone[2].x + backbone[2].radius * 0.8, backbone[2].y);
                this.ctx.stroke();
                
                this.ctx.lineWidth = 3;  // Reduced from 4
                this.ctx.beginPath();
                this.ctx.moveTo(backbone[2].x - backbone[2].radius * 0.7, backbone[2].y - 4);
                this.ctx.lineTo(backbone[3].x + backbone[3].radius * 0.7, backbone[3].y);
                this.ctx.stroke();
                
                // Draw side chain (smaller)
                this.drawOrientedSideChain(backbone[1].x, backbone[1].y, residueType);
                
                // Draw atoms
                backbone.forEach(atom => {
                    // Add subtle glow for high frequency contacts
                    if (frequency > 0.7) {
                        this.ctx.save();
                        this.ctx.shadowBlur = 10;
                        this.ctx.shadowColor = atom.color;
                        this.ctx.fillStyle = atom.color;
                        this.ctx.beginPath();
                        this.ctx.arc(atom.x, atom.y, atom.radius, 0, 2 * Math.PI);
                        this.ctx.fill();
                        this.ctx.restore();
                    }
                    
                    this.ctx.fillStyle = atom.color;
                    this.ctx.strokeStyle = 'black';
                    this.ctx.lineWidth = 1.5;  // Reduced from 2
                    
                    this.ctx.beginPath();
                    this.ctx.arc(atom.x, atom.y, atom.radius, 0, 2 * Math.PI);
                    this.ctx.fill();
                    this.ctx.stroke();
                    
                    if (atom.name === 'N' || atom.name.startsWith('O')) {
                        this.ctx.fillStyle = 'white';
                        this.ctx.font = 'bold 13px Arial';  // Reduced from 16px
                        this.ctx.textAlign = 'center';
                        this.ctx.textBaseline = 'middle';
                        const label = atom.name.startsWith('O') ? 'O' : atom.name;
                        this.ctx.fillText(label, atom.x, atom.y);
                    }
                });
                
                this.ctx.restore();
                
                // Draw residue name
                this.ctx.fillStyle = 'darkred';
                this.ctx.font = 'bold 20px Arial';  // Reduced from 26px
                this.ctx.textAlign = 'center';
                this.ctx.textBaseline = 'bottom';
                this.ctx.fillText(residueName, x, y - 38);  // Adjusted position
                
                // Draw rank badge
                this.ctx.fillStyle = '#FFD700';
                this.ctx.strokeStyle = '#000000';
                this.ctx.lineWidth = 1.5;
                this.ctx.beginPath();
                this.ctx.arc(x + 40, y - 25, 14, 0, 2 * Math.PI);  // Smaller badge
                this.ctx.fill();
                this.ctx.stroke();
                this.ctx.fillStyle = 'black';
                this.ctx.font = 'bold 14px Arial';  // Reduced from 18px
                this.ctx.textAlign = 'center';
                this.ctx.textBaseline = 'middle';
                // 计算TOP3排名 - 使用residueName参数
                const top3Contacts = this.contacts.filter(c => c.isTop3).sort((a, b) => b.frequency - a.frequency);
                const rankIndex = top3Contacts.findIndex(c => c.residue === residueName) + 1;
                this.ctx.fillText('#' + rankIndex, x + 40, y - 25);
            }
            
            drawOrientedSideChain(x, y, residueType) {
                this.ctx.strokeStyle = '#000000';
                this.ctx.lineWidth = 3;  // Reduced from 4
                
                if (residueType === 'GLY') {
                    this.ctx.fillStyle = '#808080';
                    this.ctx.font = 'bold 12px Arial';  // Reduced from 14px
                    this.ctx.textAlign = 'center';
                    this.ctx.fillText('H', x, y + 20);  // Adjusted position
                } else {
                    this.ctx.beginPath();
                    this.ctx.moveTo(x, y + 11);  // Adjusted positions
                    this.ctx.lineTo(x, y + 28);
                    this.ctx.stroke();
                    
                    this.ctx.fillStyle = '#808080';
                    this.ctx.strokeStyle = 'black';
                    this.ctx.lineWidth = 1.5;  // Reduced from 2
                    this.ctx.beginPath();
                    this.ctx.arc(x, y + 28, 8, 0, 2 * Math.PI);  // Smaller R group
                    this.ctx.fill();
                    this.ctx.stroke();
                    
                    this.ctx.fillStyle = 'white';
                    this.ctx.font = 'bold 10px Arial';  // Reduced from 12px
                    this.ctx.textAlign = 'center';
                    this.ctx.fillText('R', x, y + 28);
                }
            }
            
            drawSimpleLabel(x, y, residueType, residueName, color, index) {
                const label = residueName;
                const radius = 28;  // Slightly smaller
                
                // Add subtle shadow for depth
                this.ctx.save();
                this.ctx.shadowBlur = 5;
                this.ctx.shadowColor = 'rgba(0,0,0,0.2)';
                this.ctx.shadowOffsetX = 2;
                this.ctx.shadowOffsetY = 2;
                
                this.ctx.fillStyle = 'white';
                this.ctx.strokeStyle = color;  // Use frequency color instead of fixed color
                this.ctx.lineWidth = 3;
                this.ctx.beginPath();
                this.ctx.arc(x, y, radius, 0, 2 * Math.PI);
                this.ctx.fill();
                this.ctx.stroke();
                
                this.ctx.restore();
                
                this.ctx.fillStyle = 'darkblue';
                this.ctx.font = 'bold 13px Arial';  // Slightly smaller
                this.ctx.textAlign = 'center';
                this.ctx.textBaseline = 'middle';
                this.ctx.fillText(label, x, y);
            }
            
            draw() {
                this.ctx.save();
                this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
                
                this.ctx.translate(this.viewOffset.x, this.viewOffset.y);
                this.ctx.scale(this.zoom, this.zoom);
                
                if (this.is3DMode) {
                    this.drawGrid();
                    this.draw3DScene();
                } else {
                    this.drawGrid();
                    if (this.showConnections) this.drawConnections();
                    this.drawLigand();
                    this.drawContacts();
                }
                
                this.ctx.restore();
            }
            
            animate() {
                this.draw();
                requestAnimationFrame(() => this.animate());
            }
            
            resetPositions() {
                this.restoreInitialState();
            }
            
            rotateWheel() {
                const rotationStep = Math.PI / 6;
                const steps = 20;
                let currentStep = 0;
                
                const animate = () => {
                    if (currentStep < steps) {
                        const delta = rotationStep / steps;
                        this.wheelRotation += delta;
                        
                        this.contacts.forEach(contact => {
                            const ligandAtom = this.ligandAtoms.find(a => a.id === contact.ligandAtom);
                            const centerX = ligandAtom ? ligandAtom.x : this.centerX;
                            const centerY = ligandAtom ? ligandAtom.y : this.centerY;
                            
                            const dx = contact.x - centerX;
                            const dy = contact.y - centerY;
                            const dist = Math.sqrt(dx * dx + dy * dy);
                            const currentAngle = Math.atan2(dy, dx);
                            const newAngle = currentAngle + delta;
                            
                            contact.x = centerX + dist * Math.cos(newAngle);
                            contact.y = centerY + dist * Math.sin(newAngle);
                            contact.angle = newAngle;
                        });
                        
                        currentStep++;
                        setTimeout(animate, 20);
                    }
                };
                animate();
            }
        }
        
        let contactMap;
        window.addEventListener('DOMContentLoaded', () => { 
            contactMap = new ContactMap(); 
        });
        
        function toggle3DMode() {
            if (!contactMap.is3DMode) {
                contactMap.is3DMode = true;
                const btn = document.getElementById('mode3DBtn');
                btn.textContent = '3D Mode';
                btn.classList.add('active');
                contactMap.rotationX = 0;
                contactMap.rotationY = 0;
                contactMap.rotationZ = 0;
            } else {
                contactMap.is3DMode = false;
                const btn = document.getElementById('mode3DBtn');
                btn.textContent = '2D Mode';
                btn.classList.remove('active');
                contactMap.restoreInitialState();
            }
        }
        
        function toggleDistanceLock() {
            contactMap.distanceLocked = !contactMap.distanceLocked;
            const btn = document.getElementById('distanceLockBtn');
            if (contactMap.distanceLocked) {
                btn.textContent = '🔒 Distance Locked';
                btn.classList.remove('secondary');
            } else {
                btn.textContent = '🔓 Distance Free';
                btn.classList.add('secondary');
            }
        }
        
        function resetPositions() { contactMap.resetPositions(); }
        function rotateWheel() { contactMap.rotateWheel(); }
        function toggleConnections() { 
            contactMap.showConnections = !contactMap.showConnections; 
            const btn = document.getElementById('connectBtn');
            btn.textContent = contactMap.showConnections ? 'Hide Connections' : 'Show Connections';
        }
        function toggleHydrogens() {
            contactMap.showHydrogens = !contactMap.showHydrogens;
            const btn = document.getElementById('hydrogenBtn');
            btn.textContent = contactMap.showHydrogens ? 'Hide H atoms' : 'Show H atoms';
        }
        
        function centerView() { 
            // Calculate the bounding box of all visible elements
            let minX = Infinity, maxX = -Infinity;
            let minY = Infinity, maxY = -Infinity;
            
            // Include ligand atoms
            contactMap.ligandAtoms.forEach(atom => {
                minX = Math.min(minX, atom.x);
                maxX = Math.max(maxX, atom.x);
                minY = Math.min(minY, atom.y);
                maxY = Math.max(maxY, atom.y);
            });
            
            // Include contacts (only if visible)
            if (!contactMap.is3DMode) {
                contactMap.contacts.forEach(contact => {
                    minX = Math.min(minX, contact.x);
                    maxX = Math.max(maxX, contact.x);
                    minY = Math.min(minY, contact.y);
                    maxY = Math.max(maxY, contact.y);
                });
            }
            
            // Calculate center and size of content
            const contentWidth = maxX - minX;
            const contentHeight = maxY - minY;
            const contentCenterX = (minX + maxX) / 2;
            const contentCenterY = (minY + maxY) / 2;
            
            // Calculate zoom to fit content
            const canvasWidth = contactMap.canvas.width;
            const canvasHeight = contactMap.canvas.height;
            const padding = 50; // Padding around content
            
            const zoomX = (canvasWidth - padding * 2) / contentWidth;
            const zoomY = (canvasHeight - padding * 2) / contentHeight;
            contactMap.zoom = Math.min(zoomX, zoomY, 2); // Cap at 2x zoom
            
            // Center the view
            contactMap.viewOffset.x = canvasWidth / 2 - contentCenterX * contactMap.zoom;
            contactMap.viewOffset.y = canvasHeight / 2 - contentCenterY * contactMap.zoom;
        }
        
        function exportImage() {
            const quality = parseInt(document.getElementById('exportQuality').value);
            const filename = document.getElementById('exportFilename').value || 'contact_analysis';
            
            // Create high-resolution canvas
            const hdCanvas = document.createElement('canvas');
            hdCanvas.width = contactMap.canvas.width * quality;
            hdCanvas.height = contactMap.canvas.height * quality;
            const hdCtx = hdCanvas.getContext('2d');
            
            // Save current state
            const originalZoom = contactMap.zoom;
            const originalOffset = {...contactMap.viewOffset};
            
            // Scale up for high quality
            hdCtx.scale(quality, quality);
            
            // Temporarily replace context
            const originalCtx = contactMap.ctx;
            contactMap.ctx = hdCtx;
            
            // Enable high quality rendering
            contactMap.ctx.imageSmoothingEnabled = true;
            contactMap.ctx.imageSmoothingQuality = 'high';
            
            // Draw the scene
            contactMap.draw();
            
            // Restore original context
            contactMap.ctx = originalCtx;
            
            // Create download link
            const link = document.createElement('a');
            link.download = `${filename}_${quality}x.png`;
            link.href = hdCanvas.toDataURL('image/png', 1.0);
            link.click();
            
            // Clean up
            hdCanvas.remove();
        }
        
        function zoomIn() { 
            const newZoom = Math.min(contactMap.zoom * 1.2, 3);
            const rect = contactMap.canvas.getBoundingClientRect();
            const centerX = rect.width / 2;
            const centerY = rect.height / 2;
            
            const worldX = (centerX - contactMap.viewOffset.x) / contactMap.zoom;
            const worldY = (centerY - contactMap.viewOffset.y) / contactMap.zoom;
            
            contactMap.zoom = newZoom;
            
            contactMap.viewOffset.x = centerX - worldX * contactMap.zoom;
            contactMap.viewOffset.y = centerY - worldY * contactMap.zoom;
        }
        
        function zoomOut() { 
            const newZoom = Math.max(contactMap.zoom * 0.8, 0.3);
            const rect = contactMap.canvas.getBoundingClientRect();
            const centerX = rect.width / 2;
            const centerY = rect.height / 2;
            
            const worldX = (centerX - contactMap.viewOffset.x) / contactMap.zoom;
            const worldY = (centerY - contactMap.viewOffset.y) / contactMap.zoom;
            
            contactMap.zoom = newZoom;
            
            contactMap.viewOffset.x = centerX - worldX * contactMap.zoom;
            contactMap.viewOffset.y = centerY - worldY * contactMap.zoom;
        }
    </script>
</body>
</html>
