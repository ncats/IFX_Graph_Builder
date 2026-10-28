"""Bound memory used by the public POUNCE sheet-table JSON endpoint."""

from starlette.responses import JSONResponse


MAX_POUNCE_JSON_BYTES = 50 * 1024 * 1024


class PounceJsonSizeLimit:
    """Reject oversized JSON bodies before FastAPI builds Pydantic sheet tables."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        path = scope.get("path", "")
        root_path = scope.get("root_path", "").rstrip("/")
        if root_path and path.startswith(root_path + "/"):
            path = path[len(root_path):]
        if scope["type"] != "http" or path != "/pounce/api/validate":
            await self.app(scope, receive, send)
            return

        response = JSONResponse(
            {"detail": f"POUNCE JSON request exceeds {MAX_POUNCE_JSON_BYTES} bytes"},
            status_code=413,
        )
        headers = dict(scope.get("headers", []))
        try:
            declared_size = int(headers.get(b"content-length", b"0"))
        except ValueError:
            declared_size = 0
        if declared_size > MAX_POUNCE_JSON_BYTES:
            await response(scope, receive, send)
            return

        received_size = 0
        buffered_messages = []
        while True:
            message = await receive()
            if message["type"] == "http.request":
                received_size += len(message.get("body", b""))
                if received_size > MAX_POUNCE_JSON_BYTES:
                    await response(scope, receive, send)
                    return
                buffered_messages.append(message)
                if not message.get("more_body", False):
                    break
            elif message["type"] == "http.disconnect":
                return

        index = 0

        async def replay_receive():
            nonlocal index
            if index < len(buffered_messages):
                message = buffered_messages[index]
                index += 1
                return message
            return {"type": "http.disconnect"}

        await self.app(scope, replay_receive, send)
