function dApp_MCP = structminus(App, MCP)
dApp_MCP.my = abs(App.my - MCP.my);
dApp_MCP.Cy = abs(App.Cy - MCP.Cy);
end